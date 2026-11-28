#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")/../src/research"

find * -mindepth 1 -name '*.py' | sort | while read -r script; do
  echo
  printf "\033[1mRunning %s...\033[0m\n" "$script"
  ./"$script"
done

echo
printf "\033[1mAll scripts ran successfully.\033[0m\n"
