#!/usr/bin/env bash
# Package the current state. Usage: bash tools/pack.sh <zip-name-without-ext>
# The ZIP always contains STATUS.md, which says what is verified and what is not at this point.
set -e
cd "$(dirname "$0")/.."
name="$1"; out="/home/claude/deliver/wip/$name.zip"
rm -rf tests/out/wav tests/out/broken-missing-field.html
rm -f "$out"
zip -qr "$out" STATUS.md PROGRESS.md CHANGELOG.md REPORT.md build.mjs tools src tests shots dist orig -x "node_modules/*"
unzip -tq "$out"
echo "files: $(unzip -l "$out" | tail -1 | awk '{print $2}')  size: $(du -h "$out" | cut -f1)"
for f in lh-tl-listening.html lh-listening.html tl-listening.html; do
  a=$(unzip -p "$out" dist/$f | sha256sum | cut -c1-12); b=$(sha256sum dist/$f | cut -c1-12); [ "$a" = "$b" ] && echo "ok $f $a" || { echo "MISMATCH $f"; exit 1; }
done
