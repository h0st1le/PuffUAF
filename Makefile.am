data_files = 
etc_files = etc/volcanos.txt etc/puffrc
doc_files = license.terms INSTALL


EXTRA_DIST = $(data_files) $(etc_files) $(doc_files)

if BUILD_INSTALL_DOCS
  DOCS=
else
  DOCS=
endif

SUBDIRS = $(DOCS) man src test web 
	
install-data-local:
	$(mkinstalldirs) $(prefix)/etc
	$(INSTALL_DATA) etc/volcanos.txt $(prefix)/etc/volcanos.txt
	$(INSTALL_DATA) etc/puffrc $(prefix)/etc/puffrc

