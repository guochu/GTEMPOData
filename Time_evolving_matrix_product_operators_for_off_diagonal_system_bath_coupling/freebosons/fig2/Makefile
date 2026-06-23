all: fig.pdf

fig.pdf: fig.tex fig-inc.eps
	pdflatex fig.tex

fig.tex fig-inc.eps: fig.gnuplot
	gnuplot fig.gnuplot

.PHONY: clean
clean:
	rm -f fig.pdf fig.tex fig.aux fig.log fig-inc.eps fig-inc-eps-converted-to.pdf
