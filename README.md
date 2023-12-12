# bianca_workshop

[![Build pages](https://github.com/UPPMAX/bianca_workshop/actions/workflows/pages/pages-build-deployment/badge.svg?branch=main)](https://github.com/UPPMAX/bianca_workshop/actions/workflows/pages/pages-build-deployment)
[![Check links](https://github.com/UPPMAX/bianca_workshop/actions/workflows/check_links.yaml/badge.svg?branch=main)](https://github.com/UPPMAX/bianca_workshop/actions/workflows/check_links.yaml)
[![Check spelling](https://github.com/UPPMAX/bianca_workshop/actions/workflows/check_spelling.yaml/badge.svg?branch=main)](https://github.com/UPPMAX/bianca_workshop/actions/workflows/check_spelling.yaml)

This repository contains the source code for the [Bianca workshop](https://uppmax.github.io/bianca_workshop/)

## Credits

The website is created using
[mkdocs-material](https://squidfunk.github.io/mkdocs-material). 

## Files

Filename                           |Descriptions
-----------------------------------|------------------------------------------------------------------------------------------------------
[mlc_config.json](mlc_config.json) |Configuration of the link checker
[.wordlist.txt](.wordlist.txt)     |Whitelisted words for the spell checker, use `pyspelling -c .spellcheck.yml` to do spellcheck locally
