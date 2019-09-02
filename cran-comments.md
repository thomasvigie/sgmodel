---
title: "cran_comments.md"
output: html_document
---

## Test environments
* Windows 10 Pro
* ubuntu 14.04 (on travis-ci), R 3.5

## R CMD check results
There were no ERRORs, WARNINGs, or NOTEs.  

## Downstream dependencies
There are currently no downstream dependencies for this package.

## General comment
Running devtools::build_win() returns the following:
checking CRAN incoming feasibility ... NOTE
Maintainer: ‘Thomas Vigie <vigiethomas@gmail.com>’

Possibly mis-spelled words in DESCRIPTION:
  Hadley (14:11)
  Tauchen (13:2)
  Wickham (14:2)
  
  But they are not mis-spelled. Following the online post https://stat.ethz.ch/pipermail/r-package-devel/2016q2/000805.html
  I ignore those NOTES.
  