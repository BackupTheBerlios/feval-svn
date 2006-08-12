master.dvi: master.tex  [a-z,A-Z]*.tex
	latex master.tex 

ps: master.dvi
	dvips -o feval.ps master

pdf: master-pdf.tex  [a-z,A-Z]*.tex
	pdflatex master-pdf.tex
	mv master-pdf.pdf feval.pdf

clean: 
	rm *.aux *.log *.rel *~

all: ps pdf clean