
#!/bin/bash

taskset -c 998 nextflow run main.nf -entry $0 -c run.config -profile odap -resume -with-report -with-trace