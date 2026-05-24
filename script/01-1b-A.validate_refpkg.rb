
STDOUT.sync = true; STDERR.sync = true

require 'rake'
require 'json'
require 'digest'
require 'fileutils'

rpkg, odir, falnO, ftreO, fhmmO, ftaxO, fposO = ARGV
name = File.basename(rpkg)

### The derived files (FastTree min-evo / gamma trees, the taxtastic
### pplacer package, backbone.mfa) depend only on the refpkg source files,
### not on the query. Cache them under <rpkg>/derived/ so repeated PiPP runs
### skip the (~tens of seconds per refpkg) FastTree + taxit work.
###
### Invalidation key: MD5 of the three source files that feed the derived
### outputs (alignment + tree + hmm), stored in <rpkg>/derived/SOURCE.md5.
### If the refpkg directory is not writable, no cache is created and the
### derived files are regenerated into odir every run (original behavior).
derived  = "#{rpkg}/derived"
key_file = "#{derived}/SOURCE.md5"
src_key  = [falnO, ftreO, fhmmO].map{ |f| Digest::MD5.file(f).to_s }.join("-")

cache_valid = File.exist?(key_file) && File.read(key_file).strip == src_key

### parse LENG / NAME of fhmm (cheap; always done)
hmmlen  = "0"
hmmname = ""
open(fhmmO){ |fr|
  while l = fr.gets
    if l =~ /^LENG\s+(\d+)/
      hmmlen = $1
    elsif l =~ /^NAME\s+(\S+)/
      hmmname = $1
    end
  end
}

### output paths in odir
faln = "#{odir}/backbone.mfa"
ftre = "#{odir}/backbone.nwk"
fhmm = "#{odir}/backbone.hmm"
ftax = ftaxO ? "#{odir}/taxon.tsv" : nil
fpos = fposO ? "#{odir}/position.tsv" : nil

apdir  = "#{odir}/for_apples-2"
ftreME = "#{apdir}/backbone_min_evo.nwk"
flogME = "#{apdir}/backbone_min_evo.log"
ferrME = "#{apdir}/backbone_min_evo.err"

ppdir  = "#{odir}/for_pplacer"
ftreGM = "#{ppdir}/backbone_gamma.nwk"
flogGM = "#{ppdir}/backbone_gamma.log"
ferrGM = "#{ppdir}/backbone_gamma.err"

# {{{ def generate_derived  (the original, expensive path)
def generate_derived(odir, falnO, ftreO, fhmmO, ftaxO, fposO,
                     faln, ftre, fhmm, ftax, fpos, apdir, ppdir, ftreME, flogME, ferrME, ftreGM, flogGM, ferrGM)
  ### backbone.mfa: copy of aligned fasta with only gene ID (drop description; witch-ng workaround)
  open(faln, "w"){ |fw|
    IO.read(falnO).split(/^>/)[1..-1].each{ |ent|
      lab, *seq = ent.split("\n")
      gid = lab.split(" ")[0]
      fw.puts ">#{gid}\n#{seq*""}"
    }
  }

  ### copies of tree / hmm / taxon.tsv / position.tsv
  open(ftre, "w"){ |fw| fw.puts IO.read(ftreO) } if ftre
  open(fhmm, "w"){ |fw| fw.puts IO.read(fhmmO) } if fhmm
  open(ftax, "w"){ |fw| fw.puts IO.read(ftaxO) } if ftax
  open(fpos, "w"){ |fw| fw.puts IO.read(fposO) } if fpos

  ### minimum-evolution tree (FastTree) for APPLES-2
  mkdir_p apdir unless Dir.exist?(apdir)
  puts "### Generating minimum evolution distanced tree using FastTree for APPLES-2..."
  sh "FastTree -nosupport -nome -noml -log #{flogME} -intree #{ftre} < #{faln} > #{ftreME} 2> #{ferrME}"
  puts ""

  ### gamma tree (FastTree) + taxtastic package for pplacer (used only when input tree is IQTREE)
  mkdir_p ppdir unless Dir.exist?(ppdir)
  puts "### Generating gamma tree using FastTree for pplacer (used only when the input tree is IQTREE)..."
  sh "FastTree -nosupport -gamma -nome -mllen -log #{flogGM} -intree #{ftre} < #{faln} > #{ftreGM} 2> #{ferrGM}"
  sh "pushd #{ppdir} && taxit create -l backbone -P backbone -t #{File.basename(ftreGM)} -s #{File.basename(flogGM)} -f ../#{File.basename(faln)} --stats-type FastTree && mv backbone/* . && rmdir backbone && popd"
  puts ""
end
# }}}

if cache_valid
  ### cache hit: copy derived files into odir (skip FastTree / taxit)
  puts "### refpkg derived cache hit: #{derived} -> #{odir}"
  Dir["#{derived}/*"].each{ |path|
    base = File.basename(path)
    next if base == "SOURCE.md5"      ### key file, not a product
    next if base == "backbone.json"   ### regenerated fresh below with odir paths
    FileUtils.cp_r(path, "#{odir}/")
  }
else
  ### cache miss: generate into odir
  generate_derived(odir, falnO, ftreO, fhmmO, ftaxO, fposO,
                   faln, ftre, fhmm, ftax, fpos, apdir, ppdir, ftreME, flogME, ferrME, ftreGM, flogGM, ferrGM)

  ### populate the cache if the refpkg directory is writable; otherwise skip
  if File.writable?(rpkg)
    begin
      FileUtils.rm_rf(derived)
      FileUtils.mkdir_p(derived)
      [faln, ftre, fhmm, ftax, fpos, apdir, ppdir].compact.each{ |path|
        FileUtils.cp_r(path, "#{derived}/") if File.exist?(path)
      }
      File.write(key_file, src_key)
      puts "### cached refpkg derived files to #{derived}"
    rescue => e
      $stderr.puts "### could not write derived cache (#{e.class}: #{e.message}); continuing without cache"
    end
  else
    puts "### refpkg #{rpkg} is not writable; skipping derived cache"
  end
end

### backbone.json is always written fresh so its paths point into this run's odir
h = { name: name, refpkg: rpkg, hmmlen: hmmlen, hmmname: hmmname,
  fhmmO: fhmmO, ftaxO: ftaxO, fposO: fposO, falnO: falnO, ftreO: ftreO,
  fhmm: fhmm, ftax: ftax, fpos: fpos, faln: faln, ftre: ftre,
  ftreME: ftreME, flogME: flogME,
  ftreGM: ftreGM, flogGM: flogGM, ppdir: ppdir,
}

fjsn = "#{odir}/backbone.json"
open(fjsn, "w"){ |fw| fw.puts h.to_json }
