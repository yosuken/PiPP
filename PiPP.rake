
STDOUT.sync = true; STDERR.sync = true

# {{{ procedures
WriteBatch  = lambda do |t, jobdir, outs|
  jdir = "#{jobdir}/#{t.name.split(":")[-1]}"; mkdir_p jdir unless File.directory?(jdir)
  jnum = outs.size

  if jnum > 0
    outs.each_slice(jnum).with_index(1){ |ls, idx| ## always 1 file
      open("#{jdir}/#{t.name.split(".")[1..-1]*"."}.sh", "w"){ |fjob|
        fjob.puts ls
      }
    }
  end
end

RunBatch  = lambda do |t, jobdir, ncpu, logdir|
  jdir = "#{jobdir}/#{t.name.split(":")[-1]}"
  ldir = "#{logdir}/#{t.name.split(":")[-1]}"; mkdir_p ldir unless File.directory?(ldir)
  tdir = "#{jdir}/tmp"; mkdir_p tdir unless File.directory?(tdir)

  Dir["#{jdir}/*.sh"].sort_by{ |fin| fin.split(".")[-1].to_i }.each{ |fin| ## always 1 or 0 file
    sh "TMPDIR=#{tdir} parallel --jobs #{ncpu} --joblog #{ldir}/parallel.log <#{fin}"
  }
  open("#{ldir}/exit", "w"){ |fw| fw.puts "exit with status 0 at #{Time.now.strftime("%Y-%m-%d_%H:%M:%S")}" }
end

PrintStatus = lambda do |current, total, status, t|
  puts ""
  puts "\e[1;32m===== #{Time.now}\e[0m"
  puts "\e[1;32m===== step #{current} / #{total} (#{t.name}) -- #{status}\e[0m"
  puts ""
  # STDOUT.flush
end

CheckVersion = lambda do |commands|
  commands.each{ |command|
    str = case command
    when "ruby"
      %{which ruby && ruby -v 2>&1}
    when "hmmsearch"
      %{which hmmsearch && hmmsearch -h 2>&1 |head -n 2}
    when "mafft"
      %{which mafft && mafft -v 2>&1 |head -n 4 |tail -n 1}
    when "gappa"
      %{which gappa && gappa --version}
    when "pplacer"
      %{which pplacer && pplacer --version}
    when "parallel"
      %{which parallel && LANG=C parallel --version 2>&1 |head -n 1}
    when "witch-ng"
      %{echo #{WITCH_NG} && #{WITCH_NG} --version 2>&1}
    when "run_apples.py"
      %{which run_apples.py && run_apples.py --version 2>&1 |tail -n 1}
    when "FastTreeMP"
      %{which FastTreeMP && FastTreeMP -expert 2>&1 |head -n 1}
    when "FastTree"
      %{which FastTree && FastTree -expert 2>&1 |head -n 1}
    when "epa-ng"
      %{which epa-ng && epa-ng -v 2>&1}
    end
    puts ""
    puts "\e[1;32m===== check version: #{command}\e[0m"
    puts ""
    puts "$ #{str}"
    ### run
    puts `#{str}`
    # STDOUT.flush
  }
end
# }}} procedures


# {{{ task controller
task :default do
  # STDOUT.sync = true ### for real-time output

  ### constants
  Errmsg    = "\e[1;31mError:\e[0m"
  Warmsg    = "\e[1;35mWarning:\e[0m"
  Sccmsg    = "\e[1;32mSuccess:\e[0m"

  ### constants from arguments
  Odir         = ENV["outdir"]            ## output directory
  OdirExist    = ENV["outdir_exist"]      ## output directory exist? ("true" or "")
  Fque         = ENV["query"]             ## query (single FASTA file)
  Rpkgs        = ENV["refpkg"]            ## refpkg
  Ncpu         = ENV["ncpus"].to_i        ## num of CPUs
  Evalue       = ENV["evalue"].to_f       ## hmmsearch evalue threshold (default: 1e-5)
  EvalueDom    = ENV["evaluedom"].to_f    ## hmmsearch domain evalue threshold (default: 1e-2)

  MinSeqLen    = ENV["minseqlen"].to_i    ## minimum query length
  MinHmmLen    = ENV["minhmmlen"].to_i    ## hmmsearch minimum region hit length of hmm
  MinHmmCov    = ENV["minhmmcov"].to_f    ## minimum fraction of hmmsearch minimum region hit length of hmm per domain length
  MinAliLen    = ENV["minalilen"].to_i    ## hmmsearch minimum region hit length of ali
  MinAliCov    = ENV["minalicov"].to_f    ## minimum fraction of hmmsearch minimum region hit length of ali per gene length
  MinHmmLenDom = ENV["minhmmlendom"].to_i ## hmmsearch minimum domain hit length of hmm
  MinHmmCovDom = ENV["minhmmcovdom"].to_f ## minimum fraction of hmmsearch minimum domain hit length of hmm per domain length
  MinAliLenDom = ENV["minalilendom"].to_i ## hmmsearch minimum domain hit length of ali
  MinAliCovDom = ENV["minalicovdom"].to_f ## minimum fraction of hmmsearch minimum domain hit length of ali per gene length

  C_size       = ENV["chunk_size"].to_i   ## chunk size
  Aligner      = ENV["aligner"]           ## witch-ng|mafft-add
  M_mafft      = ENV["mafft_method"]      ## E-INS-i|FFT-NS-i|FFT-NS-2|FFT-NS-i_addfragments|FFT-NS-2_addfragments
  EPA_NG_model = ENV["epa_ng_model"]      ## required if --placer epa-ng is given
  Ex_lvs       = ENV["extract_levels"]    ## clade/taxonomy level for 'gappa prepare extract' (default: 0)
  Placer       = ENV["placer"]            ## pplacer|apples-2|epa-ng
  KeepIntermediate = ENV["keep_intermediate"] == "true"  ## keep prefilter/chunks/batch/log/tasks after run
  Z            = 1_000_000                ## hmmsearch data size for evalue calculation
  # Trim_opt     = ENV["trim_option"]   ## 'merge' (merge regions of hmmserach hits)  or 'largest' (take largest hit) (default: merge)

  ### define tasks
  tasks = []
  tasks << "01-1a-A.validate_query"
  tasks << "01-1a-B.parse_query_info"
  tasks << "01-1b-A.validate_refpkg"
  tasks << "01-1b-B.parse_refpkg_info"
  tasks << "01-2a.hmmsearch"
  tasks << "01-2b.parse_hmmsearch"
  tasks << "01-2c.split_fasta"
  tasks << "01-2d.copy_detected"

  if ENV["only_detect"] != "true"
    tasks << "01-2e.prepare_for_placement"
    tasks << "01-3a.chunkify"

    case Aligner
    when "witch-ng"
      tasks << "01-3b.witch-ng"
    when "mafft-add"
      tasks << "01-3b.mafft-add"
    else
      raise("#{Errmsg} --aligner #{Aligner} is not available. Choose either 'witch-ng' or 'mafft-add'.")
    end

    case Placer
    when "pplacer"
      tasks << "01-3c.pplacer"
    when "apples-2"
      tasks << "01-3c.apples-2"
    when "epa-ng"
      tasks << "01-3c.epa-ng"
    else
      raise("#{Errmsg} --placer #{Placer} is not available. Choose either 'pplacer', 'apples-2' or 'epa-ng'.")
    end

    tasks << "01-3d.unchunkify"
    tasks << "01-3e.unchunkify_alignment"
    tasks << "01-4a.info"
    tasks << "01-4b.lwr-list"
    tasks << "01-4d.assign"  ## if ftax is given
    tasks << "01-4f.graft"
    tasks << "01-4h.aligned_position" ## if fpos is given
    tasks << "01-6a.aa_feature"
    tasks << "01-7a.duckdb_import"

    # tasks << "01-4e.extract" ## if ftax is given ### memory requirement is too high for rep_2022-06-08
    # tasks << "01-4g.heat-tree"
    # tasks << "01-5a.krd"

    ### currently not used
    ### tasks << "01-4c.edpl"    ## high memory requirement. --> do not run
    ### tasks << "01-5b.edgepca" ## high memory requirement if input (num of samples?) is large.
    ### tasks << "01-5c.squash"
    ### tasks << "01-5d.dispersion"
  end

  ### store Tasks
  Tasks = tasks

  ### check witch-ng
  WITCH_NG = "#{__dir__}/bin/witch-ng"

  ### locate pipp_util (bundled Rust binary).
  ### Search order: env PIPP_UTIL_BIN > $PATH > <project>/rust/target/release/pipp_util
  PIPP_UTIL = (
    if (e = ENV["PIPP_UTIL_BIN"]) and File.executable?(e)
      e
    elsif !(p = `command -v pipp_util 2>/dev/null`.chomp).empty?
      p
    elsif File.executable?(d = "#{__dir__}/rust/target/release/pipp_util")
      d
    else
      raise "#{Errmsg} pipp_util (bundled Rust binary) not found.\n" \
            "         Install via Bioconda, or build from source:\n" \
            "         cargo build --release --manifest-path #{__dir__}/rust/Cargo.toml"
    end
  )

  ### check version
  commands  = %w|hmmsearch mafft gappa pplacer ruby|
  commands += %w|parallel witch-ng FastTree run_apples.py epa-ng|
  CheckVersion.call(commands)

  ## dir path
  Jobdir    = "#{Odir}/batch"
  Pkgdir    = "#{Odir}/refpkg"              # "#{Odir}/prefilter/refpkgs/refpkg_info.tsv
  Predir    = "#{Odir}/prefilter" 
  PreQuedir = "#{Odir}/prefilter/query"     # "#{Odir}/prefilter/query/#{que[:name]}/*.fa" (use .fa extension)
  PreFildir = "#{Odir}/prefilter/hmmsearch" # "#{Odir}/prefilter/hmmsearch/#{que[:name]}/{evalues.tsv,*.domtableout}"
  Cnkdir    = "#{Odir}/chunks"              # "#{Odir}/chunks/#{pkg[:name]}/{chunk,alignment,placement}"
  Resdir    = "#{Odir}/result"              # "#{Odir}/result/#{pkg[:name]}/{seq,alignment,placement,assign,extract,graft,feature}"
  Logdir    = "#{Odir}/log/tasks"

  ## variables
  $fque     = nil
  $fpkgs    = []
  $rnames   = {}
  $ex_lvs   = []

  ## validate extract_levels
  ex_lvs_errmsg = "#{Errmsg} --extract-levels '#{Ex_lvs}' should be either -1, 0, positive integer, or comma separeted positive integers."
  case Ex_lvs
  when "-1", /^[0-9,]+$/
    a = Ex_lvs.split(",")
    raise(ex_lvs_errmsg) if a.size > 1 and a.include?("0")
    $ex_lvs = a.map(&:to_i).uniq.sort
  else raise(ex_lvs_errmsg)
  end

  ## Odir exist?
  $stderr.puts "\n\n#{Warmsg} output directory #{Odir} already exists. Overwrite it.\n\n" if OdirExist == "true"

  ### run
  NumStep  = Tasks.size
  Tasks.each.with_index(1){ |task, idx|
    Rake::Task[task].invoke(idx)
    # STDOUT.flush
  }

  ### cleanup intermediate dirs (suppress with --keep-intermediate)
  unless KeepIntermediate
    puts "\n\e[1;32m===== cleanup intermediate files\e[0m"
    [Predir, Cnkdir, Jobdir, Logdir].each{ |d| rm_rf d if File.directory?(d) }
  end
end
# }}} default (run all tasks)


# {{{ tasks
# {{{ desc "01-1a-A.validate_query"
desc "01-1a-A.validate_query"
task "01-1a-A.validate_query", ["step"] do |t, args|
  ### [!!!] should not be skipped even if already done
  PrintStatus.call(args.step, NumStep, "START", t)
  outs   = []

  ## validate query and make a copy
  fque = Fque
  raise("#{Errmsg} <query> file: #{fque} is empty.") if File.zero?(fque)

  if fque =~ /\.gz$/ or fque =~ /\.gzip/
    require 'zlib'
    Zlib::GzipReader.open(fque) { |fr|
      l = fr.gets
      raise("#{Errmsg} <query> file: #{fque} has wrong format.\nThe first line should begin '>'") if l[0] != ">"
    }
  else
    File.open(fque){ |fr|
      l = fr.gets
      raise("#{Errmsg} <query> file: #{fque} has wrong format.\nThe first line should begin '>'") if l[0] != ">"
    }
  end

  $stderr.puts ["", "", "\e[1;32m===== check query file \e[0m"]

  ### parse fasta and check sequence length
  mkdir_p PreQuedir
  name = File.basename(fque).gsub(/\.gz/, "").gsub(/\.gzip/, "").split(".")[0..-2]*"."
  fa   = "#{PreQuedir}/#{name}.fa"
  fjsn = "#{PreQuedir}/#{name}.json"
  flst = "#{PreQuedir}/#{name}.list"

  if !File.exist?(fjsn)
    outs << "#{PIPP_UTIL} validate-query --min-seq-len #{MinSeqLen} --name #{name} -i #{fque} --fa #{fa} --json #{fjsn} --list #{flst}"
  end

  next if outs.size == 0

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
# }}}

# {{{ desc "01-1a-B.parse_query_info"
desc "01-1a-B.parse_query_info"
task "01-1a-B.parse_query_info", ["step"] do |t, args|
  ### [!!!] should not be skipped even if already done
  PrintStatus.call(args.step, NumStep, "START", t)
  require 'json'

  name = File.basename(Fque).gsub(/\.gz/, "").gsub(/\.gzip/, "").split(".")[0..-2]*"."
  fjsn = "#{PreQuedir}/#{name}.json"
  sjsn = IO.readlines(fjsn)[0]
  puts sjsn ### print info

  $fque = JSON.parse(sjsn, symbolize_names: true)
end
# }}}

# {{{ desc "01-1b-A.validate_refpkg"
desc "01-1b-A.validate_refpkg"
task "01-1b-A.validate_refpkg", ["step"] do |t, args|
  ### [!!!] should not be skipped even if already done
  PrintStatus.call(args.step, NumStep, "START", t)
  outs   = []
  script = "#{__dir__}/script/#{t.name}.rb"

  ## valdiate refpkg (does not genrate local copy file)
  ### [2022-06-15] sort files
  rpkgs = Rpkgs.split(/[,\s]+/).sort_by{ |path| File.basename(path) }.inject([]){ |a, path| a += Dir[path.gsub("~", ENV["HOME"])].sort }
  $stderr.puts ["", "", "\e[1;32m===== check refpkg (N=#{rpkgs.size}) \e[0m"]
  raise("#{Errmsg} no refpkg directory detected.") if rpkgs.size == 0

  ### refpkg name duplication check
  $rnames = {} ### refpkg names
  rpkgs.each{ |rpkg|
    name = File.basename(rpkg)
    raise("#{Errmsg} refpkg name #{name} is not unique.") if $rnames[name]
    $rnames[name] = 1

    fa    = Dir["#{rpkg}/*.fa"] + Dir["#{rpkg}/*.mfa"] + Dir["#{rpkg}/*.fasta"] + Dir["#{rpkg}/*.faa"]
    fhmm  = Dir["#{rpkg}/*.hmm"]
    ftre  = Dir["#{rpkg}/*.tree"] + Dir["#{rpkg}/*.nwk"] + Dir["#{rpkg}/*.newick"]
    ftax  = Dir["#{rpkg}/taxon.tsv"] 
    fpos  = Dir["#{rpkg}/position.tsv"]

    raise("#{Errmsg} #{rpkg} is not a directory.") unless File.directory?(rpkg)
    raise("#{Errmsg} #{rpkg} does not contain fasta file. #{rpkg}/*{.fa|.mfa|.faa|.fasta} should exist.") if fa.size == 0
    raise("#{Errmsg} #{rpkg} contains multiple fasta files. #{rpkg}/*{.fa|.mfa|.faa|.fasta} should be only one.") if fa.size > 1
    raise("#{Errmsg} #{rpkg} contains multiple tree files. #{rpkg}/*{.tree|.nwk|.newick} should be only one.") if ftre.size > 1
    raise("#{Errmsg} #{rpkg} does not contain .hmm file. #{rpkg}/*.hmm should exist.") if fhmm.size == 0
    raise("#{Errmsg} #{rpkg} contains multiple .hmm files. #{rpkg}/*.hmm should be only one.") if fhmm.size > 1
    raise("#{Errmsg} #{fhmm[0]} contains .hmm file with multiple hmms. HMM should be only one per a file.") if IO.read(fhmm[0]).split(/\/\/\s+/).size > 1

    ### [!!!] TODO: sequence (not name) in tree should be nonredundant
    ### [!!!] TODO: sequence name in tree should be same as sequence name in alignment
    ### [!!!] TODO: ftax format check
    ### [!!!] TODO: fpos format check

    falnO  = fa[0]
    ftreO  = ftre[0]
    fhmmO  = fhmm[0]
    ftaxO  = ftax.size == 0 ? nil : ftax[0] ### optional
    fposO  = fpos.size == 0 ? nil : fpos[0] ### optional

    odir = "#{Pkgdir}/#{name}"; mkdir_p odir unless File.directory?(odir)
    faln = "#{Pkgdir}/#{name}/backbone.mfa"
    flog = "#{Pkgdir}/#{name}/backbone.log"
    if !File.exist?(faln)
      outs << "ruby #{script} #{rpkg} #{odir} #{falnO} #{ftreO} #{fhmmO} #{ftaxO} #{fposO} >#{flog} 2>&1"
    end
  }

  next if outs.size == 0

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
# }}}

# {{{ desc "01-1b-B.parse_refpkg_info"
desc "01-1b-B.parse_refpkg_info"
task "01-1b-B.parse_refpkg_info", ["step"] do |t, args|
  ### [!!!] should not be skipped even if already done
  PrintStatus.call(args.step, NumStep, "START", t)
  require 'json'

  $rnames.each_key{ |name|
    fjsn = "#{Pkgdir}/#{name}/backbone.json"
    sjsn = IO.readlines(fjsn)[0]
    puts sjsn ### print info

    # $fpkgs << eval(sjsn)
    $fpkgs << JSON.parse(sjsn, symbolize_names: true)
  }
end
# }}}

# {{{ desc "01-2a.hmmsearch"
desc "01-2a.hmmsearch"
task "01-2a.hmmsearch", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit") ### skip if already done
  outs   = []

  npara = [Ncpu, $fpkgs.size].min
  ncpu  = [1, (Ncpu.to_f / npara).round].max

  puts "###"
  puts "### number of parallel for hmmsearch: #{npara}"
  puts "### number of cpu per a hmmsearch: #{ncpu} CPUs (total: #{Ncpu})"
  puts "###"

  p $fque
  raise("#{Errmsg} query #{$fque[:name]} has no sequence.") if $fque[:numseq] == 0
  odir = "#{PreFildir}/#{$fque[:name]}/out"; mkdir_p odir
  $fpkgs.each{ |pkg|
    db   = pkg[:fhmm]
    fa   = $fque[:fasta]
    fout = "#{odir}/#{pkg[:name]}.out"
    flog = "#{odir}/#{pkg[:name]}.log"
    outs << "hmmsearch --cpu #{ncpu} -Z #{Z} --notextw -o #{fout} #{db} #{fa} >#{flog} 2>&1"
  }

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, npara, Logdir)
end
# }}}

# {{{ desc "01-2b.parse_hmmsearch"
desc "01-2b.parse_hmmsearch"
task "01-2b.parse_hmmsearch", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit") ### skip if already done
  outs   = []

  fa    = $fque[:fasta]
  idir  = "#{PreFildir}/#{$fque[:name]}/out"

  if Dir["#{idir}/*.out"].size > 0
    ### concat hmmsearch output files
    odir = "#{PreFildir}/#{$fque[:name]}/parsed"; mkdir_p odir unless File.directory?(odir)
    fhmm = "#{odir}/hmmsearch_concat.out"
    open(fhmm, "w"){ |fw|
      $fpkgs.each{ |pkg|
        f = "#{idir}/#{pkg[:name]}.out"
        if File.exist?(f) and !File.zero?(f)
          fw.puts IO.read(f)
        else
          raise("#{Errmsg} #{f} does not exist or zero-sized.")
        end
      }
    }

    ### [2026-05-23 use pipp_util parse-hmmsearch]
    flog    = "#{odir}/parse.log"
    option  = "--gene-evalue #{Evalue} --evalue #{EvalueDom} --create-evalue-table --min-hmm-len-dom #{MinHmmLenDom} --min-hmm-cov-dom #{MinHmmCovDom} --min-ali-len-dom #{MinAliLenDom} --min-ali-cov-dom #{MinAliCovDom} "
    option += "--min-hmm-len #{MinHmmLen} --min-hmm-cov #{MinHmmCov} --min-ali-len #{MinAliLen} --min-ali-cov #{MinAliCov}"
    outs << "#{PIPP_UTIL} parse-hmmsearch #{option} -i #{fhmm} -f #{fa} -o #{odir} >#{flog} 2>&1"
  end

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
# }}}

# {{{ desc "01-2c.split_fasta"
desc "01-2c.split_fasta"
task "01-2c.split_fasta", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit") ### skip if already done

  # {{{ def parse_fasta(fin, fout, ids)
  def parse_fasta(fin, fout, ids)
    open(fout, "w"){ |fw|
      flag, ent = 0, []
      open(fin){ |fr|
        while l = fr.gets
          if l[0] == ">"
            fw.puts ent if ent != [] 
            flag, ent = 0, []

            id = l.strip[1..-1].split(/\s+/)[0]
            if ids[id]
              flag = 1
              ent << l.strip
            end
          else
            ent << l.strip if flag == 1
          end
        end
        fw.puts ent if ent != []
      }
    }
  end
  # }}}

  ### parse parse_hmmsearch.rb output
  ftsv = "#{PreFildir}/#{$fque[:name]}/parsed/best-hit.tsv"
  hmm2regs = {} ### hmm_name => regions
  hmm2gids = {} ### hmm_name => gene_ids
  open(ftsv){ |fr|
    _ = fr.gets # skip header
    #       0           1             2         3        4         5        6      7     8         9
    # protein  length(aa)  protein_info  hmm_name  hmm_acc  hmm_desc  hmm_len  score  bias  c-Evalue
    #       10      11      12      13      14      15      16   17           18          19
    # i-Evalue  hmm.fm  hmm.to  ali.fm  ali.to  env.fm  env.to  acc  full-Evalue  full-score
    #   20           21
    # link  region_name
    while l = fr.gets
      a = l.strip.split("\t", -1)
      gid, hmm, reg = a.values_at(0, 3, 21)
      hmm2gids[hmm] ||= {}
      hmm2gids[hmm][gid] = 1
      hmm2regs[hmm] ||= {}
      hmm2regs[hmm][reg] = 1
    end
  }

  $fpkgs.each{ |pkg|
    regs = hmm2regs[pkg[:hmmname]]
    gids = hmm2gids[pkg[:hmmname]]
    next if !regs or regs.keys.size == 0

    %w|whole region|.zip(%w|best-hit.whole.fa best-hit.fa|, [gids, regs]){ |type, fname, ids|
      fin  = "#{PreFildir}/#{$fque[:name]}/parsed/#{fname}"
      odir = "#{PreFildir}/#{$fque[:name]}/seq/#{type}/#{pkg[:name]}"; mkdir_p odir unless File.directory?(odir)
      fout = "#{odir}/#{$fque[:name]}.fa"

      parse_fasta(fin, fout, ids)
      # p [fin, fout, ids.size]
    }
  }
end
# }}}

# {{{ desc "01-2d.copy_detected"
desc "01-2d.copy_detected"
task "01-2d.copy_detected", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit") ### skip if already done
  outs = []

  $fpkgs.each{ |pkg|
    %w|whole region|.each{ |type|
      fin  = "#{PreFildir}/#{$fque[:name]}/seq/#{type}/#{pkg[:name]}/#{$fque[:name]}.fa"
      next unless File.exist?(fin)

      odir = "#{Resdir}/#{pkg[:name]}/seq"; mkdir_p odir unless File.directory?(odir)
      fout = "#{odir}/#{type}.fa"
      outs << "cp #{fin} #{fout}"
    }
  }

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
# }}}

# {{{ desc "01-2e.prepare_for_placement"
desc "01-2e.prepare_for_placement"
task "01-2e.prepare_for_placement", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)

  ### assign Npara and NcpuP by counting number of refpkgs with queries identified.
  npkg = 0
  $fpkgs.each{ |pkg|
    fas = "#{PreFildir}/#{$fque[:name]}/seq/region/#{pkg[:name]}/#{$fque[:name]}.fa"
    next if Dir[fas].size == 0
    npkg += 1
  }

  Npara = [Ncpu, npkg].min

  ### number of cpu per a refpkg
  NcpuP = [1, (Ncpu.to_f / npkg).round].max

  puts "###"
  puts "### number of parallel for placement: #{npkg}"
  puts "### number of cpu per a placement: #{NcpuP} CPUs (total: #{Ncpu})"
  puts "###"
end
# }}}

# {{{ desc "01-3a.chunkify"
desc "01-3a.chunkify"
task "01-3a.chunkify", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit") ### skip if already done

  outs = []
  cmd  = "gappa prepare chunkify --threads 1 --allow-file-overwriting"

  $fpkgs.each{ |pkg|
    fas = "#{PreFildir}/#{$fque[:name]}/seq/region/#{pkg[:name]}/#{$fque[:name]}.fa"
    next if Dir[fas].size == 0

    ### prepare input fasta file (cat all query files)
    odir = "#{Cnkdir}/#{pkg[:name]}/chunk"; mkdir_p odir
    flog = "#{odir}/chunkify.log"
    outs << "#{cmd} --chunk-size #{C_size} --chunks-out-dir #{odir} --abundances-out-dir #{odir} --fasta-path #{fas} >#{flog} 2>&1"
  }

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
# }}}

# {{{ desc "01-3b.witch-ng"
desc "01-3b.witch-ng"
task "01-3b.witch-ng", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit") ### skip if already done
  outs = []

  HmmSizeLb = 100 ### default: 10 (it may take long time.)

  script0 = "#{__dir__}/script/extract_ori_reg.rb" ## extract region of the input alignment, since witch-ng output is sometimes longer than input
  script1 = "#{__dir__}/script/UO2X.rb" ## "U" --> "X" in aligned fasta
  script2 = "#{__dir__}/script/only_query.rb" ## "U" --> "X" in aligned fasta

  $fpkgs.each{ |pkg|
    fas = Dir["#{Cnkdir}/#{pkg[:name]}/chunk/chunk_*.fasta"].sort_by{ |i| File.basename(i).gsub(/^chunk_/, "").gsub(/\.fasta$/, "").to_i }
    next if fas.size == 0

    ### [!!!] faln should be nonredundant (sequences are same as tree) --> TODO: validation process

    fas.each{ |fa|
      chnk  = File.basename(fa).gsub(/\.fasta$/, "") ## chunk_0, chunk_1, ...
      odir  = "#{Cnkdir}/#{pkg[:name]}/alignment/#{chnk}"; mkdir_p odir
      flog  = "#{odir}/witch-ng.log"
      ftmp0 = "#{odir}/witch-ng.orig.fa" ## original output (sometimes longer than input)
      fout  = "#{odir}/witch-ng.fa"
      ftmp1 = "#{odir}/witch-ng.UO2X.fa" ## used only by pplacer ('U' --> 'X')
      ftmp2 = "#{odir}/witch-ng.UO2X.only_query.fa"

      # witch-ng add --threads 4 -i queries.fa -b backbone.afa -t backbone.tre -o extended_alignment.afa
      out   = []
      out  << "RUST_BACKTRACE=full #{WITCH_NG} add --threads #{NcpuP} --hmm-size-lb #{HmmSizeLb} -i #{fa} -b #{pkg[:faln]} -t #{pkg[:ftre]} -o #{ftmp0} >#{flog} 2>&1" ## compatible with "U"
      out  << "ruby #{script0} #{pkg[:faln]} #{ftmp0} #{fout}" ### extract region of the input alignment, since witch-ng output is sometimes longer than input
      out  << "ruby #{script1} #{fout} #{ftmp1}"
      out  << "ruby #{script2} #{pkg[:faln]} #{ftmp1} #{ftmp2}"
      outs << out*" && "
    }
  }

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Npara, Logdir)
end
# }}}

# {{{ desc "01-3b.mafft-add"
desc "01-3b.mafft-add"
task "01-3b.mafft-add", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit") ### skip if already done
  outs = []

  script1 = "#{__dir__}/script/UO2X.rb" ## "U" --> "X" in aligned fasta
  script2 = "#{__dir__}/script/only_query.rb" ## "U" --> "X" in aligned fasta

  $fpkgs.each{ |pkg|
    fas = Dir["#{Cnkdir}/#{pkg[:name]}/chunk/chunk_*.fasta"].sort_by{ |i| File.basename(i).gsub(/^chunk_/, "").gsub(/\.fasta$/, "").to_i }
    next if fas.size == 0

    ### [!!!] faln should be nonredundant (sequences are same as tree) --> TODO: validation process

    ### make clean header for faln
    faln = "#{Cnkdir}/#{pkg[:name]}/backbone.mfa"
    open(faln, "w"){ |fw|
      IO.read(pkg[:faln]).split(/^>/)[1..-1].each{ |ent|
        lab, *seq = ent.split("\n")
        gid = lab.split(" ")[0]
        fw.puts ">#{gid}\n#{seq*""}"
      }
    }

    fas.each{ |fa|
      chnk  = File.basename(fa).gsub(/\.fasta$/, "") ## chunk_0, chunk_1, ...
      odir  = "#{Cnkdir}/#{pkg[:name]}/alignment/#{chnk}"; mkdir_p odir
      flog  = "#{odir}/mafft-add.log"
      fout  = "#{odir}/mafft-add.fa"
      ftmp1 = "#{odir}/mafft-add.UO2X.fa" ## used only by pplacer ('U' --> 'X')
      ftmp2 = "#{odir}/witch-ng.UO2X.only_query.fa"

      ## [!!!] --maxiterate > 2 is not compatible with --keeplength
      if M_mafft == "E-INS-i"
        option = "--genafpair" ## slow and accurate
        add    = "--add"
      elsif M_mafft == "FFT-NS-i" 
        option = "--maxiterate 2"  ## faster
        add    = "--add"
      elsif M_mafft == "FFT-NS-2" 
        option = "--retree 2"  ## fastest
        add    = "--add"
      elsif M_mafft == "FFT-NS-i_addfragments" 
        option = "--maxiterate 2"  ## faster
        add    = "--addfragments"
      elsif M_mafft == "FFT-NS-2_addfragments" 
        option = "--retree 2"  ## fastest
        add    = "--addfragments"
      else
        raise("#{Errmsg} --mafft-method #{M_mafft} is not available.")
      end

      out   = []
      out  << "mafft #{option} --anysymbol --thread #{NcpuP} #{add} #{fa} --keeplength #{faln} >#{fout} 2>#{flog}" ## compatible with "U"
      out  << "ruby #{script1} #{fout} #{ftmp1}"
      out  << "ruby #{script2} #{pkg[:faln]} #{ftmp1} #{ftmp2}"
      outs << out*" && "

      ## [?] it could be automatic selection if option == ""
    }
  }

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Npara, Logdir)
end
# }}}

# {{{ desc "01-3c.pplacer"
desc "01-3c.pplacer"
task "01-3c.pplacer", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit") ### skip if already done
  outs = []

  $fpkgs.each{ |pkg|
    ### Aligner: mafft-add or witch-ng
    fas = Dir["#{Cnkdir}/#{pkg[:name]}/alignment/chunk_*/#{Aligner}.UO2X.fa"] ### ref and query fasta
    fas = fas.sort_by{ |fa| fa.split("/")[-2].gsub(/^chunk_/, "").to_i }

    is_valid_pkg_given = true ### false if (1) CONTENTS.json does not exist, or (2) pplacer with IQ-TREE tree is given

    fas.each{ |fa|
      ### mafft-add.UO2X.fa or witch-ng.UO2X.fa
      chnk  = fa.split("/")[-2]
      fcont = "#{pkg[:refpkg]}/CONTENTS.json"

      if File.exist?(fcont)
        # parse a line like below
        # "phylo_model": "phylo_modelxmpfgc8r.json",
        require 'json'
        h = JSON.parse(IO.read(fcont))
        raise("#{Errmsg} #{fcont} does not contain phylo_model info.") unless h["files"] and h["files"]["phylo_model"]
        fphylo = "#{pkg[:refpkg]}/#{h["files"]["phylo_model"]}"

        # parse a line like below
        # "program": "IQ-TREE 2.2.2.6",
        h = JSON.parse(IO.read(fphylo))
        raise("#{Errmsg} #{fphylo} does not contain program info.") unless h["program"]
        prog = h["program"].split(" ")[0]

        if prog == "IQ-TREE"
          is_valid_pkg_given = false
          puts "#{Warmsg} pplacer with IQ-TREE tree is not supported. Run pplacer alternatively with FastTree gamma distance tree for #{pkg[:name]}."
        end
      else
        is_valid_pkg_given = false
        puts "#{Warmsg} #{fcont} does not exist. Run pplacer with FastTree gamma distance tree for #{pkg[:name]}."
      end

      odir  = "#{Cnkdir}/#{pkg[:name]}/placement/#{chnk}"; mkdir_p odir
      flog  = "#{odir}/pplacer.log"

      if is_valid_pkg_given
        outs << "pplacer --verbosity 2 -j #{NcpuP} --out-dir #{odir} -c #{pkg[:refpkg]} #{fa} >#{flog} 2>&1"
      else
        outs << "pplacer --verbosity 2 -j #{NcpuP} --out-dir #{odir} -c #{pkg[:ppdir]} #{fa} >#{flog} 2>&1"
      end
    }
  }

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Npara, Logdir)
end
# }}}

# {{{ desc "01-3c.apples-2"
desc "01-3c.apples-2"
task "01-3c.apples-2", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit") ### skip if already done
  outs = []

  $fpkgs.each{ |pkg|
    ### Aligner: mafft-add or witch-ng
    fas = Dir["#{Cnkdir}/#{pkg[:name]}/alignment/chunk_*/#{Aligner}.UO2X.fa"] ### ref and query fasta
    fas = fas.sort_by{ |fa| fa.split("/")[-2].gsub(/^chunk_/, "").to_i }

    fas.each{ |fa|
      ### mafft-add.UO2X.fa or witch-ng.UO2X.fa
      chnk  = fa.split("/")[-2]
      odir  = "#{Cnkdir}/#{pkg[:name]}/placement/#{chnk}"; mkdir_p odir
      flog  = "#{odir}/apples-2.log"
      fout  = "#{odir}/#{File.basename(fa).gsub(/\.fa$/, "")}.jplace"

      outs << "run_apples.py -T #{NcpuP} -o #{fout} -s #{pkg[:faln]} -t #{pkg[:ftreME]} -x #{fa} >#{flog} 2>&1"
    }
  }

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Npara, Logdir)
end
# }}}

# {{{ desc "01-3c.epa-ng"
desc "01-3c.epa-ng"
task "01-3c.epa-ng", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit") ### skip if already done
  outs = []

  $fpkgs.each{ |pkg|
    ### Aligner: mafft-add or witch-ng
    fas = Dir["#{Cnkdir}/#{pkg[:name]}/alignment/chunk_*/#{Aligner}.UO2X.only_query.fa"]  ### only query fasta
    fas = fas.sort_by{ |fa| fa.split("/")[-2].gsub(/^chunk_/, "").to_i }

    fas.each{ |fa|
      ### mafft-add.UO2X.fa or witch-ng.UO2X.fa
      chnk  = fa.split("/")[-2]
      fpkg  = pkg[:refpkg] ### [!!!] should be nonredundant (sequences are same as tree)
      odir  = "#{Cnkdir}/#{pkg[:name]}/placement/#{chnk}"; mkdir_p odir
      flog  = "#{odir}/epa-ng.log"
      _fout = "#{odir}/epa_result.jplace"
      fout  = "#{odir}/#{Aligner}.UO2X.jplace"

      # epa-ng --ref-msa $REF_MSA --tree $TREE --query $QRY_MSA --model $MODEL
      # EPA_NG_model = "LG" 
      # EPA_NG_model = "PROTGTR" 
      outs << "epa-ng -T #{NcpuP} -m #{EPA_NG_model} -w #{odir} -s #{pkg[:faln]} -t #{pkg[:ftre]} -q #{fa} >#{flog} 2>&1 && mv #{_fout} #{fout}"
    }
  }

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Npara, Logdir)
end
# }}}

# {{{ desc "01-3d.unchunkify"
desc "01-3d.unchunkify"
task "01-3d.unchunkify", ["step"] do |t, args|
  ### [!!!] [2020-11-17] Use of parallel cause "core dump". DO NOT USE parallel
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit") ### skip if already done
  outs = []
  cmd1 = "gappa prepare unchunkify --threads 1 --allow-file-overwriting"

  $fpkgs.each{ |pkg|
    fplcs = "#{Cnkdir}/#{pkg[:name]}/placement/chunk_*/*.jplace"
    next if Dir[fplcs].size == 0

    fabus = "#{Cnkdir}/#{pkg[:name]}/chunk/abundances_*.json"
    odir  = "#{Resdir}/#{pkg[:name]}/placement"; mkdir_p odir
    flog  = "#{odir}/unchunkify.log"

    outs << "#{cmd1} --out-dir #{odir} --abundances-path #{fabus} --jplace-path #{fplcs} >#{flog} 2>&1"
  }

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
# }}} 

# {{{ desc "01-3e.unchunkify_alignment"
desc "01-3e.unchunkify_alignment"
task "01-3e.unchunkify_alignment", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit") ### skip if already done
  outs = []

  script = "#{__dir__}/script/#{t.name}.rb"

  $fpkgs.each{ |pkg|
    fas  = "#{Cnkdir}/#{pkg[:name]}/chunk/chunk_*.fasta" ## region seq, hashed
    alns = "#{Cnkdir}/#{pkg[:name]}/alignment/chunk_*/#{Aligner}.fa" ## aligned seq, hashed

    ques = Dir["#{PreFildir}/#{$fque[:name]}/seq/region/#{pkg[:name]}/#{$fque[:name]}.fa"] ## seq_region, labeled

    next if ques.size == 0

    adir = "#{Resdir}/#{pkg[:name]}/alignment"; mkdir_p adir unless File.directory?(adir)

    outs << "ruby #{script} #{pkg[:faln]} '#{fas}' '#{alns}' '#{ques*","}' #{adir}"
  }

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
# }}}

# {{{ desc "01-4a.info"
desc "01-4a.info"
task "01-4a.info", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit") ### skip if already done
  outs = []
  cmd  = "gappa examine info --threads 1 --allow-file-overwriting"

  $fpkgs.each{ |pkg|
    bdir  = "#{Resdir}/#{pkg[:name]}"
    fplcs = Dir["#{bdir}/placement/*.jplace"][0]
    next unless fplcs

    flog  = "#{bdir}/placement/info.log"
    outs << "#{cmd} --jplace-path #{fplcs} >#{flog} 2>&1"
  }

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
# }}}

# {{{ desc "01-4b.lwr-list"
desc "01-4b.lwr-list"
task "01-4b.lwr-list", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit") ### skip if already done
  outs = []
  cmd  = "gappa examine lwr-list --threads 1 --allow-file-overwriting"

  $fpkgs.each{ |pkg|
    bdir  = "#{Resdir}/#{pkg[:name]}"
    fplcs = Dir["#{bdir}/placement/*.jplace"][0]
    next unless fplcs

    odir  = "#{bdir}/placement"
    flog  = "#{odir}/lwr-list.log"
    outs << "#{cmd} --out-dir #{odir} --jplace-path #{fplcs} >#{flog} 2>&1"
  }

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
# }}}

# {{{ desc "01-4c.edpl"
desc "01-4c.edpl"
task "01-4c.edpl", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit") ### skip if already done
  outs = []
  cmd  = "gappa examine edpl --no-list-file --threads 1 --allow-file-overwriting"
  ### [!!!] --no-list-file: need to limit memory usage.
  ### If set, do not write out the EDPL per pquery, but just the histogram file. As the list needs to keep all pquery names in memory (to get the correct order), the memory requirements might be too large. In that case, this option can help.

  ## [!!!] memory requirement is very high even when --no-list-file is used.
  $stderr.puts "Currently, do not execute 'gappa examine edpl' due to high memory requirement."
  next ## SKIP THIS TASK !!!

  $fpkgs.each{ |pkg|
    bdir  = "#{Resdir}/#{pkg[:name]}"
    fplcs = Dir["#{bdir}/placement/*.jplace"][0]
    next unless fplcs

    odir  = "#{bdir}/placement"
    flog  = "#{odir}/edpl.log"
    pref  = "edpl_"
    outs << "#{cmd} --out-dir #{odir} --file-prefix #{pref} --jplace-path #{fplcs} >#{flog} 2>&1"
  }

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
# }}}

# {{{ desc "01-4d.assign"
desc "01-4d.assign"
task "01-4d.assign", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  # (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit") ### skip if already done

  outs = []
  # cmd  = "gappa examine assign --threads 1 --per-query-results --krona --allow-file-overwriting" ### for v0.7.1
  cmd  = "gappa examine assign --threads 1 --krona --allow-file-overwriting" ### for v0.6.0

  $fpkgs.each{ |pkg|
    bdir = "#{Resdir}/#{pkg[:name]}"
    ftax = pkg[:ftax]

    ($stderr.puts "For #{pkg}, taxon.tsv is not found. Skip the taxonomy/clade assignment step."; next) unless ftax

    fplcs = Dir["#{bdir}/placement/*.jplace"][0]
    next unless fplcs

    odir  = "#{bdir}/assign"; mkdir_p odir
    flog  = "#{odir}/assign.log"
    outs << "#{cmd} --out-dir #{odir} --jplace-path #{fplcs} --taxon-file #{ftax} --per-query-results >#{flog} 2>&1"
  }

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
# }}}

# {{{ desc "01-4e.extract"
desc "01-4e.extract"
task "01-4e.extract", ["step"] do |t, args|
  next if $ex_lvs == [-1] ### do not extract

  PrintStatus.call(args.step, NumStep, "START", t)
  # (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit") ### skip if already done

  outs = []
  cmd  = "gappa prepare extract --threads 1 --allow-file-overwriting"

  $fpkgs.each{ |pkg|
    bdir = "#{Resdir}/#{pkg[:name]}"
    ftax = pkg[:ftax]

    ($stderr.puts "For #{pkg}, taxon.tsv is not found. Skip the taxonomy/clade assignment step."; next) unless ftax

    ### parse levels
    max_lv = IO.readlines(ftax).map{ |l| i = l.chomp.split("\t")[1]; i ? i.split(/\s*;\s*/).size : 0 }.sort[-1]
    next if max_lv == 0
    lvs = if $ex_lvs == [0] then (1..max_lv).to_a
    else $ex_lvs.select{ |i| i <= max_lv }
    end

    fplcs = Dir["#{bdir}/placement/*.jplace"][0]
    next unless fplcs

    lvs.each{ |lv|
      fas   = "#{bdir}/seq/whole.fa"
      odir  = "#{bdir}/extract/level#{lv}"; mkdir_p odir
      flog  = "#{odir}/extract.log"

      fcld  = "#{odir}/clade.tsv"
      open(fcld, "w"){ |fw|
        fw.puts IO.readlines(ftax).map{ |l| a = l.chomp.split("\t"); [a[0], a[1].split(/\s*;\s*/)[lv - 1]]*"\t" }
      }

      outs << "#{cmd} --samples-out-dir #{odir} --sequences-out-dir #{odir} --jplace-path #{fplcs} --fasta-path #{fas} --clade-list-file #{fcld} >#{flog} 2>&1"
    }
  }

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
# }}}

# {{{ desc "01-4f.graft"
desc "01-4f.graft"
task "01-4f.graft", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit") ### skip if already done
  outs = []
  cmd  = "gappa examine graft --threads 1 --fully-resolve --name-prefix Q_ --allow-file-overwriting"

  $fpkgs.each{ |pkg|
    bdir  = "#{Resdir}/#{pkg[:name]}"
    fplcs = Dir["#{bdir}/placement/*.jplace"][0]
    next unless fplcs

    odir  = "#{bdir}/graft"; mkdir_p odir
    flog  = "#{odir}/graft.log"
    outs << "#{cmd} --out-dir #{odir} --jplace-path #{fplcs} >#{flog} 2>&1"
  }

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
# }}}

# {{{ desc "01-4g.heat-tree"
desc "01-4g.heat-tree"
task "01-4g.heat-tree", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit") ### skip if already done
  outs = []
  opt_svg = "--svg-tree-shape circular --color-list viridis --reverse-color-list --svg-tree-stroke-width 3 --svg-tree-ladderize"
  prefix  = "tree"
  cmd  = "gappa examine heat-tree --threads 1 #{opt_svg} --tree-file-prefix #{prefix} --write-newick-tree --write-nexus-tree --write-phyloxml-tree --write-svg-tree --allow-file-overwriting"

  script = "#{__dir__}/script/nexus2itol.rb"

  $fpkgs.each{ |pkg|
    bdir  = "#{Resdir}/#{pkg[:name]}"
    fplcs = Dir["#{bdir}/placement/*.jplace"][0]
    next unless fplcs

    odir  = "#{bdir}/heat-tree"; mkdir_p odir
    flog  = "#{odir}/heat-tree.log"
    outs << "#{cmd} --out-dir #{odir} --jplace-path #{fplcs} >#{flog} 2>&1 && ruby #{script} #{odir}/#{prefix}.nexus #{odir}/#{prefix}.itol.txt"
  }

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
# }}}

# {{{ desc "01-4h.aligned_position"
desc "01-4h.aligned_position"
task "01-4h.aligned_position", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  # (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit") ### skip if already done
  outs = []

  script = "#{__dir__}/script/#{t.name}.rb"

  $fpkgs.each{ |pkg|
    bdir    = "#{Resdir}/#{pkg[:name]}"
    fpos    = pkg[:fpos]
    frefaln = pkg[:faln]

    ($stderr.puts "For #{pkg}, position.tsv is not found. Skip to check alignment positions."; next) unless fpos

    fplcs = Dir["#{bdir}/placement/*.jplace"][0]
    next unless fplcs

    odir    = "#{bdir}/alignment"
    fallaln = "#{odir}/aligned.fa" ## should exist.
    fasn    = "#{bdir}/assign/per_query.tsv" ## use result of assignment. If not exist, do not use assignment result
    flog    = "#{odir}/aligned_position.log"
    outs << "ruby #{script} #{odir} #{fpos} #{fallaln} #{frefaln} #{fasn} >#{flog} 2>&1"
  }

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
# }}}

# {{{ desc "01-5a.krd"
desc "01-5a.krd"
task "01-5a.krd", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit") ### skip if already done
  outs = []
  cmd  = "gappa analyze krd --threads 1 --allow-file-overwriting"

  $fpkgs.each{ |pkg|
    bdir = "#{Resdir}/#{pkg[:name]}"

    fplcs = "#{bdir}/placement/*.jplace"
    next if Dir[fplcs].size < 2

    odir  = "#{bdir}/krd"; mkdir_p odir
    flog  = "#{odir}/krd.log"
    outs << "#{cmd} --out-dir #{odir} --jplace-path #{fplcs} >#{flog} 2>&1"
  }

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
# }}}

# {{{ desc "01-5b.edgepca"
desc "01-5b.edgepca"
task "01-5b.edgepca", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit") ### skip if already done
  outs = []
  cmd  = "gappa analyze edgepca --threads 1 --write-nexus-tree --write-phyloxml-tree --write-svg-tree --allow-file-overwriting"
  ## newick tree is same as tree in refpkg (no color information) --> do not use --write-newick-tree

  $fpkgs.each{ |pkg|
    bdir = "#{Resdir}/#{pkg[:name]}"

    fplcs = "#{bdir}/placement/*.jplace"
    next if Dir[fplcs].size < 2

    odir  = "#{bdir}/edgepca"; mkdir_p odir
    flog  = "#{odir}/edgepca.log"
    outs << "#{cmd} --out-dir #{odir} --jplace-path #{fplcs} >#{flog} 2>&1"
  }

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
# }}}

# {{{ desc "01-5c.squash"
desc "01-5c.squash"
task "01-5c.squash", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit") ### skip if already done
  outs = []
  cmd  = "gappa analyze squash --threads 1 --write-nexus-tree --write-phyloxml-tree --write-svg-tree --allow-file-overwriting"
  ## newick tree is same as tree in refpkg (no color information) --> do not use --write-newick-tree
  ## [!!!] this task will generate (2n - 2) * 3 tree files where n is num of input query fasta files.

  $fpkgs.each{ |pkg|
    bdir = "#{Resdir}/#{pkg[:name]}"

    fplcs = "#{bdir}/placement/*.jplace"
    next if Dir[fplcs].size < 2

    odir  = "#{bdir}/squash"; mkdir_p odir
    flog  = "#{odir}/squash.log"
    outs << "#{cmd} --out-dir #{odir} --jplace-path #{fplcs} >#{flog} 2>&1"
  }

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
# }}}

# {{{ desc "01-5d.dispersion"
desc "01-5d.dispersion"
task "01-5d.dispersion", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit") ### skip if already done
  outs = []
  cmd  = "gappa analyze dispersion --threads 1 --write-nexus-tree --write-phyloxml-tree --write-svg-tree --allow-file-overwriting"
  ## newick tree is same as tree in refpkg (no color information) --> do not use --write-newick-tree

  $fpkgs.each{ |pkg|
    bdir = "#{Resdir}/#{pkg[:name]}"

    fplcs = "#{bdir}/placement/*.jplace"
    next if Dir[fplcs].size < 2

    odir  = "#{bdir}/dispersion"; mkdir_p odir
    flog  = "#{odir}/dispersion.log"
    outs << "#{cmd} --out-dir #{odir} --mass-norm absolute --jplace-path #{fplcs} >#{flog} 2>&1"
  }

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
# }}}

# {{{ desc "01-6a.aa_feature"
desc "01-6a.aa_feature"
task "01-6a.aa_feature", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit") ### skip if already done
  outs = []

  script = "#{__dir__}/script/#{t.name}.rb"

  $fpkgs.each{ |pkg|
    bdir = "#{Resdir}/#{pkg[:name]}"
    fa = "#{PreFildir}/#{$fque[:name]}/seq/region/#{pkg[:name]}/#{$fque[:name]}.fa"
    next unless File.exist?(fa)

    odir  = "#{bdir}/feature/aa"; mkdir_p odir unless File.directory?(odir)
    fout  = "#{odir}/feature.tsv"
    flog  = "#{odir}/feature.log"
    outs << "ruby #{script} #{fout} #{fa} >#{flog} 2>&1"
  }

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
# }}}

# {{{ desc "01-7a.duckdb_import"
desc "01-7a.duckdb_import"
task "01-7a.duckdb_import", ["step"] do |t, args|
  PrintStatus.call(args.step, NumStep, "START", t)
  (puts "Already done. Skipped." ; next) if File.exist?("#{Logdir}/#{t.name.split(":")[-1]}/exit") ### skip if already done
  outs = []

  $fpkgs.each{ |pkg|
    bdir = "#{Resdir}/#{pkg[:name]}"
    next unless File.directory?(bdir)

    flog = "#{bdir}/pipp.duckdb.log"
    outs << "#{PIPP_UTIL} import #{bdir} --refpkg #{pkg[:name]} --overwrite >#{flog} 2>&1"
  }

  WriteBatch.call(t, Jobdir, outs)
  RunBatch.call(t, Jobdir, Ncpu, Logdir)
end
# }}}
# }}} tasks


# require "minitest/test_task"
# Minitest::TestTask.create # named test, sensible defaults
# Minitest::TestTask.new do |t|
#   t.libs << "test"
#   t.libs << "lib"
#   t.libs << "."
#   t.warning = false
#   t.test_globs = ["test/**/*_test.rb"]
# end
# task :default => :test
