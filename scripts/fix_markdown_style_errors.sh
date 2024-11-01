#!/bin/bash
#
# Fix markdown style errors,
# as is recommended by the Markdown checker.
#
# Usage:
#
#   ./scripts/fix_markdown_style_errors.sh

if [[ "$PWD" =~ scripts$ ]]; then
    echo "FATAL ERROR."
    echo "Please run the script from the project root. "
    echo "Present working director: $PWD"
    echo " "
    echo "Tip: like this"
    echo " "
    echo "  ./scripts/fix_markdown_style_errors.sh"
    echo " "
    exit 42
fi

markdownlint --fix "**/lesson_plans/*.md"
markdownlint --fix "**/reflections/*.md"
markdownlint --fix "*.md"
markdownlint --fix "docs/intermediate/transfer.md"
markdownlint --fix "docs/intermediate/summary.md"
markdownlint --fix "docs/intermediate/ides.md"
markdownlint --fix "docs/intermediate/overview2.md"
markdownlint --fix "docs/intermediate/intro.md"
markdownlint --fix "docs/intermediate/complex_jobs.md"
markdownlint --fix "docs/intermediate/efficient_jobs.md"
markdownlint --fix "docs/intermediate/replicate_jobs.md"
markdownlint --fix "docs/intermediate/evaluation_intermediate.md"
# markdownlint --fix "**/docs/*.md"




