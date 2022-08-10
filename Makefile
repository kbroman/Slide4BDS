R_OPTS=--no-save --no-restore --no-init-file --no-site-file

broman_slide.pdf: broman_slide.tex Figs/lod_curves.png
	xelatex $<

Figs/lod_curves.png: R/lod_curves.R
	cd R;R $(R_OPTS) -e "source('$(<F)', echo=TRUE)"
