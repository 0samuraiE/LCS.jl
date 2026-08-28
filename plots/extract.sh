#!/bin/sh

set -eu

usage() {
    echo "Usage: $0 <font-size> [destination-directory]" >&2
    exit 2
}

[ "$#" -ge 1 ] && [ "$#" -le 2 ] || usage

font_size=$1
case "$font_size" in
    ''|*[!0-9]*)
        echo "font-size must be a positive integer: $font_size" >&2
        exit 2
        ;;
esac
[ "$font_size" -gt 0 ] || {
    echo "font-size must be a positive integer: $font_size" >&2
    exit 2
}

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
repo_dir=$(dirname -- "$script_dir")
source_dir="$repo_dir/figs"
destination_dir=${2:-"$repo_dir/figs-$font_size"}

mkdir -p -- "$destination_dir"

copied=0
skipped=0
for figure in \
    "$source_dir"/*-"$font_size".pdf \
    "$source_dir"/*-strong-"$font_size"-1.pdf \
    "$source_dir"/*-weak-"$font_size"-1.pdf
do
    [ -f "$figure" ] || continue

    destination="$destination_dir/$(basename -- "$figure")"
    if [ -e "$destination" ]; then
        echo "skip: $destination already exists"
        skipped=$((skipped + 1))
        continue
    fi

    cp -- "$figure" "$destination"
    copied=$((copied + 1))
done

echo "copied $copied figure(s) for font size $font_size to $destination_dir; skipped $skipped"
