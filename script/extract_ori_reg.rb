
# Extract, from a witch-ng output alignment, only the columns that correspond
# to the reference (backbone) alignment — witch-ng sometimes adds query
# insertion columns, making its output wider than the backbone.
#
# Columns are identified by content: a backbone column is the string of
# residues across the backbone sequences; the matching output column has the
# identical string. The backbone sequences are located in the output BY NAME
# (not by position), because witch-ng may emit them in a different order
# (e.g. when reusing a prebuilt eHMM via `-b <ehmm-dir>`).

faln, fin, fout = ARGV

N_check = 10000 ### compare using at most this many backbone sequences

# {{{ def read_fasta(path)  -> [ordered_names, name=>seq] (single or multi-line)
def read_fasta(path)
  order = []
  seq   = {}
  cur   = nil
  open(path){ |fr|
    while l = fr.gets
      if l[0] == ">"
        cur = l.strip[1..-1].split(/\s+/)[0]
        order << cur
        seq[cur] = +""
      elsif cur
        seq[cur] << l.strip
      end
    end
  }
  [order, seq]
end
# }}}

### backbone alignment
bb_order, bb_seq = read_fasta(faln)
bb_order = bb_order.first(N_check)
N = bb_order.size
L = bb_seq[bb_order[0]].size

### witch-ng output (indexed by name so sequence order does not matter)
_out_order, out_seq = read_fasta(fin)

### the same backbone sequences, taken from the output by name, in backbone order
out_bb = bb_order.map{ |nm|
  out_seq[nm] or raise("Error: backbone sequence #{nm} not found in witch-ng output #{fin}")
}
N1 = out_bb.size
L1 = out_bb[0].size
raise("Error: unexpected alignment. number of backbone sequences in input alignment (#{N}) and witch-ng output (#{N1}) should be the same") if N != N1

### columns as strings (over the backbone sequences, backbone order)
aln0 = (0...L).map{ |i| bb_order.map{ |nm| bb_seq[nm][i] }.join }
aln1 = (0...L1).map{ |i| out_bb.map{ |s| s[i] }.join }

### map each output column to a backbone column (monotonic; backbone columns
### appear in order, separated by query-insertion columns)
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

### write every output sequence, keeping only the backbone-corresponding columns
open(fout, "w"){ |fw|
  open(fin){ |fr|
    while l = fr.gets
      if l[0] == ">"
        fw.puts l
      else
        seq = l.strip.split(//)
        fw.puts seq.values_at(*pos1).join("")
      end
    end
  }
}
