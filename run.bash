#!/bin/bash
grep . $@ \
  |  grep ":#BSUB \|CPU time\|Max Memory\|Average Memory\|Delta Memory\|Max Swap\|Max Processes\|Max Threads\|Run time\|Turnaround time" \
  |  sed 's/:/@/' \
  |  tr -s ' ' '@' \
  |  sed 's/@/_/2' \
  |  replace_with_tab '@' \
  |  sed 's/@:@/@/' \
  |  tr -d '#' \
  |  sed 's/@/#/' \
  |  sed 's/Average_Memory/Mem_Ave/' \
  |  sed 's/Max_Memory/Mem_Max/' \
  |  sed 's/Delta_Memory/Mem_Delta/' \
  |  sed 's/\/.command.log//g' \
  |  sed 's/nf-//' \
  |  sed 's/work\///' \
  |  grep -v "#BSUB_-o" \
  |  sort \
  |  mergeuniq -merge  \
  |  sed 's/"//g' \
  |  sed 's/@rusage/ BSUB_-R#rusage/' \
  |  tr -d '[>=])' \
  |  sed 's/selectmem//' \
  |  sed 's/rusagemem//' \
  |  sed 's/spanhosts//' \
  |  tr ' #' '\t' \
  |  recut 3,1,5,7,11,13,15,17,19,23,25,27,29,31,33,35,21,37,9 \
  |  tr '@' ' ' \
  |  sed 's/_(/\t/' \
  |  sed 's/sec./sec/g' \
  |  perl -lane 'BEGIN{print "process\ttag\thash\tmemory\tcpus\trselect\trusage\tspan\tsla\ttime\tprocesses\tswap\tthreads\tmem_ave\tmem_delta\tmem_max\trun_time\tcpu_time\ttotal_time\tworkdir"};print $_'
