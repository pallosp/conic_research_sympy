#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")/../src/research"

for script in */*.py; do
  echo
  printf "\033[1mRunning $script...\033[0m"
  echo
  ./"$script"
done

echo
printf "\033[1mAll scripts ran successfully.\033[0m\n"
