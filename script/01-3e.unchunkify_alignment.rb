
STDOUT.sync = true; STDERR.sync = true

fori, fas, alns, ques, adir = ARGV

## fori: alignment given by refpkg
### fas  = "#{Cnkdir}/#{pkg[:name]}/chunk/chunk_*.fasta"
fas  =  Dir[fas].sort_by{ |fa| fa.split("/")[-1].split(".")[0].gsub(/^chunk_/, "").to_i }
### alns = "#{Cnkdir}/#{pkg[:name]}/alignment/chunk_*/#{Aligner}.fa"
alns = Dir[alns].sort_by{ |fa| fa.split("/")[-2].gsub(/^chunk_/, "").to_i }
### ques = "#{PreFildir}/#{$fque[:name]}/seq/region/#{pkg[:name]}/#{$fque[:name]}.fa"
ques = ques.split(",")

hash2seq = {}
fas.each{ |fa| ## region seq, hashed
  IO.read(fa).split(/^>/)[1..-1].each{ |ent|
    hash, *seq = ent.split("\n")
    seq = seq.join.gsub(/\s+/, "")
    hash2seq[hash] = seq
  }
}

seq2aln = {}
alns.each{ |fa| ## alingned seq, hashed
  IO.read(fa).split(/^>/)[1..-1].each{ |ent|
    hash, *aln = ent.split("\n")
    aln = aln.join.gsub(/\s+/, "")
    seq = hash2seq[hash]
    seq2aln[seq]  = aln
  }
}

header = [] ### alignment given by refpkg
IO.read(fori).split(/^>/)[1..-1].each{ |ent|
  lab, *seq = ent.split("\n")
  seq = seq.join.gsub(/\s+/, "")
  header << ">#{lab}" << seq
}

fwa0 = open("#{adir}/aligned.fa", "w")        ## alignment with ref
fwa1 = open("#{adir}/aligned_wo_ref.fa", "w") ## alignment without ref
fwa0.puts header

ques.each{ |fa| ## region seq, labeled
  IO.read(fa).split(/^>/)[1..-1].each{ |ent|
    lab, *seq = ent.split("\n")
    seq = seq.join.gsub(/\s+/, "")
    aln = seq2aln[seq]

    fwa0.puts [">"+lab, aln]
    fwa1.puts [">"+lab, aln]
  }
}

[fwa0, fwa1].each{ |fwa| fwa.close }
