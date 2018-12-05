#!/usr/bin/env bash

rsync -avz --partial \
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
           --exclude 'nohup*' \
           --exclude 'graphs/*.log' \
           evan@palaeoprime:~/qpbrute/ ~/Dropbox/Code/qpbrute/
