#!/usr/bin/env bash
# Pre-commit hook: ensure docs/api.md is up to date with pydoc-markdown.

CONFIG_FILE="pydoc-markdown.yml"
API_FILE="docs/api.md"

if [[ ! -f "$CONFIG_FILE" ]]; then
  echo "Config file $CONFIG_FILE not found" >&2
  exit 1
fi

if [[ ! -f "$API_FILE" ]]; then
  echo "$API_FILE does not exist. Run 'pydoc-markdown' to create it." >&2
  exit 1
fi

# Create a temp config to write to stdout
# It must be in the project root directory for the module discover to work.
TMP_CONFIG=".tmp.pydoc-markdown.yml"
trap 'rm -f "$TMP_CONFIG"' EXIT
grep -v '^[[:space:]]*filename:' "$CONFIG_FILE" > "$TMP_CONFIG"

# Diff generated docs against committed version
if ! .venv/bin/pydoc-markdown "$TMP_CONFIG" | cmp -s "$API_FILE" -; then
  echo
  echo "$API_FILE is out of date. Run 'pydoc-markdown' to regenerate." >&2
  exit 1
fi
