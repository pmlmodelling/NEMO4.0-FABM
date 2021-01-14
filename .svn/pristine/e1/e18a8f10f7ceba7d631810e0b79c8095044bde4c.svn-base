**********************
NEMO Quick Start Guide
**********************

| The NEMO guide is made up of several files written in
  `ReStructuredText <http://docutils.sourceforge.net/rst.html>`_ (`.rst` extension),
  a WYSIWYG markup language used in the Python community, and scattered all over the NEMO sources.
| You can view them one by one in plain text from `./source` folder, or export all to a user-friendly guide under `./build` (only HTML format at the moment, PDF expected later).

Build and export the guide in HTML
==================================

1.  Install Sphinx documentation generator, its BibTeX extension and "Read The Docs" theme thanks to `pip` packages tool

    .. code-block:: console

        $ pip install sphinx sphinxcontrib.bibtex sphinx_rtd_theme

2.  Build the HTML export with `make` in `./build/html`

    .. code-block:: console

        $ make html

3.  Finally browse the guide by opening `./build/html/NEMO_guide.html`


Edit the sources and check the output in real time
==================================================

| To facilitate the update of the guide, editors can install a useful package that will automatically trigger a new build and the reload of the HTML page for every recorded change in the sources.
| So the reviewer saves time by controlling on-line their modifications almost as it types and also by avoiding repeated interactive rebuilds.

Install `sphinx-autobuild` package

.. code-block:: console

	$ pip install sphinx-autobuild

Launch a local web server hosting a draft export of the guide (build this time in `./build/livehtml`)

.. code-block:: console

    $ make livehtml

| Open in the same time the 2 formats of the content to review: the source file and the web page by browsing from the new guide hosted by your local server on `<http://127.0.0.1:8000/NEMO_guide.html>`_.
| Start the update, save your changes and verify instantly the HTML export in your browser.

.. warning::

    | Your modifications are not taken into account?
    | For symlink file, you will have to close it to update the HTML export. Otherwise look at the log of the Sphinx build, you probably made a typo!

.. hint::

    Are there broken links? Fix "Page not found" errors by running `make linkcheck`
