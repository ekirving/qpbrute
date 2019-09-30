#!/usr/bin/env bash

rsync -az  --partial \
           --progress \
           --stats \
           --exclude '.git*' \
           --exclude '*.sh' \
           --exclude '*.txt' \
           --exclude '*.py*' \
           --exclude '*.r' \
           --exclude '*.R' \
           --exclude '*.old*' \
           --exclude '*tmp*' \
           --exclude '*.gz' \
           --exclude '*.bak*' \
           --exclude 'nohup*' \
           --exclude 'graphs*/*' \
           --exclude 'bayes*/*' \
           evan@palaeoprime:~/qpbrute/ ~/Dropbox/Code/qpbrute/
