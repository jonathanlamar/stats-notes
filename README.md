# Stats Notes

I am taking notes on a stats class on edX [here](https://courses.edx.org/courses/course-v1:MITx+18.6501x+3T2019/course/).  Since I have sworn never to compile a LaTeX document again, I am taking notes in markdown and compiling with pandoc.  Compile with the command
```bash
./mdpdf -i <FILE>.md -o <FILE>.pdf
```
If you want to run a script in the background that recompiles every time you save, run
```bash
./pandoc_loop -f <FILE WITHOUT EXTENSION>
```
These should work in Linux as long as you install the underlying utilities in those scripts.
