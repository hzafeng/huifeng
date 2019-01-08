pandoc -V mainfont="SimSun" -f markdown -t rst ./source/Ecoli.md -o ./source/Ecoli.rst
make html
git add * -f
git commit -m first
git push first master
