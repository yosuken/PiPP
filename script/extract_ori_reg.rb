
# Post-process a witch-ng output alignment:
#   (1) keep only the columns that correspond to the reference (backbone)
#       alignment — witch-ng adds query insertion columns, making its output
#       wider than the backbone;
#   (2) emit the sequences in a deterministic, canonical order — backbone
#       sequences in backbone.mfa order, then query sequences in the chunk
#       input order. witch-ng reorders its output (especially when reusing a
#       prebuilt eHMM via `-b <ehmm-dir>`), so without this the official
#       output order would depend on the eHMM cache; reordering makes it
#       identical regardless.
#
# usage: ruby extract_ori_reg.rb <backbone.mfa> <chunk.fasta> <witch-ng.out> <out.fa>

faln, fque, fin, fout = ARGV

N_check = 10000 ### compare using at most this many backbone sequences

# {{{ def read_fasta(path) -> [ordered_names, name=>seq, name=>header_after_gt]
def read_fasta(path)
  order = []
  seq   = {}
  hdr   = {}
  cur   = nil
  open(path){ |fr|
    while l = fr.gets
      if l[0] == ">"
        h   = l.strip[1..-1]
        cur = h.split(/\s+/)[0]
        order << cur
        seq[cur] = +""
        hdr[cur] = h
      elsif cur
        seq[cur] << l.strip
      end
    end
  }
  [order, seq, hdr]
end
# }}}

### backbone alignment
bb_order, bb_seq, = read_fasta(faln)
bb_order = bb_order.first(N_check)
N = bb_order.size
L = bb_seq[bb_order[0]].size

### query order (chunk input order)
q_order, = read_fasta(fque)

### witch-ng output (indexed by name; its own order is ignored)
_out_order, out_seq, out_hdr = read_fasta(fin)

### same backbone sequences taken from the output by name, in backbone order
out_bb = bb_order.map{ |nm|
  out_seq[nm] or raise("Error: backbone sequence #{nm} not found in witch-ng output #{fin}")
}
N1 = out_bb.size
L1 = out_bb[0].size
raise("Error: unexpected alignment. number of backbone sequences in input alignment (#{N}) and witch-ng output (#{N1}) should be the same") if N != N1

### columns as strings (over the backbone sequences, backbone order)
aln0 = (0...L).map{ |i| bb_order.map{ |nm| bb_seq[nm][i] }.join }
aln1 = (0...L1).map{ |i| out_bb.map{ |s| s[i] }.join }

### map each output column to a backbone column (monotonic)
conv = {} ### aln1 (output) position -> aln0 (backbone) position
k = 0
(0..L1-1).each{ |i|
  (k..L-1).each{ |j|
    if aln1[i] == aln0[j]
      conv[i] = j
      k = j + 1
      break
    end
  }
}

if conv.size != L
  raise("Error: only #{conv.size} positions could be mapped out of #{L}, when parsing #{faln} and #{fin}")
end

pos1 = conv.keys ### output column indices that correspond to backbone columns (ascending)

### write in canonical order: backbone (backbone order), then queries (chunk order),
### keeping only the backbone-corresponding columns.
open(fout, "w"){ |fw|
  (bb_order + q_order).each{ |nm|
    s = out_seq[nm] or raise("Error: sequence #{nm} not found in witch-ng output #{fin}")
    fw.puts ">#{out_hdr[nm]}"
    fw.puts pos1.map{ |i| s[i] }.join("")
  }
}
