#!/bin/bash
set -euo pipefail

cd "$(dirname "$0")"

for script in *.py; do
  [ "$script" = "__init__.py" ] && continue
  echo
  echo -e "\033[1mRunning $script...\033[0m"
  echo
  ./"$script"
done

echo
echo "\033[1mAll scripts ran successfully.\033[0m"
