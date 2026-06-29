
STDOUT.sync = true; STDERR.sync = true

require 'rake'
require 'json'
require 'digest'
require 'fileutils'

rpkg, odir, falnO, ftreO, fhmmO, ftaxO, fposO, only_detect = ARGV
falnO = nil if falnO == "__nil__"
ftreO = nil if ftreO == "__nil__"
fhmmO = nil if fhmmO == "__nil__"
ftaxO = nil if ftaxO == "__nil__"
fposO = nil if fposO == "__nil__"
only_detect = only_detect == "true"
name = File.basename(rpkg)

raise "Error: missing HMM file argument" unless fhmmO && File.exist?(fhmmO)

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

if only_detect
  FileUtils.mkdir_p(odir)
  fhmm = "#{odir}/backbone.hmm"
  ftax = ftaxO ? "#{odir}/taxon.tsv" : nil
  fpos = fposO ? "#{odir}/position.tsv" : nil

  FileUtils.cp(fhmmO, fhmm)
  FileUtils.cp(ftaxO, ftax) if ftaxO
  FileUtils.cp(fposO, fpos) if fposO

  h = { name: name, refpkg: rpkg, hmmlen: hmmlen, hmmname: hmmname,
    fhmmO: fhmmO, ftaxO: ftaxO, fposO: fposO, falnO: falnO, ftreO: ftreO,
    fhmm: fhmm, ftax: ftax, fpos: fpos, faln: falnO, ftre: ftreO,
    ftreME: nil, flogME: nil,
    ftreGM: nil, flogGM: nil, ppdir: nil,
  }

  fjsn = "#{odir}/backbone.json"
  open(fjsn, "w"){ |fw| fw.puts h.to_json }
  exit
end

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
### Bump DERIV_VERSION whenever generate_derived's logic changes so existing
### caches are invalidated even though the source files are unchanged.
###   v2: clamp negative branch lengths in the min-evo tree (apples-2 NaN fix)
DERIV_VERSION = "2"
src_key  = ([DERIV_VERSION] + [falnO, ftreO, fhmmO].map{ |f| Digest::MD5.file(f).to_s }).join("-")

cache_valid = File.exist?(key_file) && File.read(key_file).strip == src_key

### Where the derived files for this run live (and what backbone.json
### references). By default, when the cache is valid we reference it *in
### place* under <rpkg>/derived/ — no per-run copy (the eHMM caches alone are
### tens of MB per refpkg, so copying them into every run dir is pure waste;
### downstream tasks only read these files). Pass --copy-refpkg (ENV
### copy_refpkg) to materialize a full isolated copy under odir instead.
### On a cache miss we must generate into a writable place: odir.
copy_refpkg = ENV["copy_refpkg"] == "true"
dest = (cache_valid && !copy_refpkg) ? derived : odir

### derived-file paths (under dest)
faln = "#{dest}/backbone.mfa"
ftre = "#{dest}/backbone.nwk"
fhmm = "#{dest}/backbone.hmm"
ftax = ftaxO ? "#{dest}/taxon.tsv" : nil
fpos = fposO ? "#{dest}/position.tsv" : nil

apdir  = "#{dest}/for_apples-2"
ftreME = "#{apdir}/backbone_min_evo.nwk"
flogME = "#{apdir}/backbone_min_evo.log"
ferrME = "#{apdir}/backbone_min_evo.err"

ppdir  = "#{dest}/for_pplacer"
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
  ### Clamp negative branch lengths (least-squares artifacts) to 0. On trees
  ### with negative branches APPLES-2 can emit NaN pendant lengths, producing
  ### a jplace that is invalid JSON and unparseable by gappa. (The output
  ### jplace tree is additionally clamped by `pipp_util clamp-jplace`.)
  me_tree = File.read(ftreME)
  File.write(ftreME, me_tree.gsub(/:-[0-9.eE+-]+/, ":0"))
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
  if dest == odir
    ### --copy-refpkg: materialize a full isolated copy under odir
    puts "### refpkg derived cache hit: copying #{derived} -> #{odir} (--copy-refpkg)"
    Dir["#{derived}/*"].each{ |path|
      base = File.basename(path)
      next if base == "SOURCE.md5"      ### key file, not a product
      next if base == "backbone.json"   ### regenerated fresh below with odir paths
      FileUtils.cp_r(path, "#{odir}/")
    }
  else
    ### default: reference the cache in place (no copy)
    puts "### refpkg derived cache hit: referencing #{derived} in place (pass --copy-refpkg to copy into the run dir)"
  end
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
