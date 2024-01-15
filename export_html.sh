#!/bin/bash

prac=$1

jupyter nbconvert --to html ceg2707_practical${prac}.ipynb --HTMLExporter.theme=dark
