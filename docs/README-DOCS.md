# How to Write Documentation

The documentation uses the [MyST](https://mystmd.org/) markup language and renders to [Sphinx](https://www.sphinx-doc.org) as HTML, PDF or EPUB documents.

In order to build the documentation locally, from your BeamTracking.jl/docs directory, 
create a software environment using the commands
```{code} bash
conda env create -y -f conda.yml    # This only needs to be done once.
conda activate myst                 # Do this with any new window.
```

Compile HTML or LaTex with the commands:
```{code} bash
make html
make latexpdf
```

Open the file `build/html/index.html` with your web browser to visualize.
You are now ready to edit the markdown files that compose the documentation!

If you like to atuomatically rebuild changes on save of edited files, 
run in your BeamTracking.jl/docs directory:
```{code} bash
sphinx-autobuild src build/html
```
and open the URL shown in the terminal, usually [http://127.0.0.1:8000](http://127.0.0.1:8000).

If you add new markdown files, do not forget to add them to the table of contents defined `docs/src/index.md`.
Finally, once you are happy with your changes, do not forget to commit and push them to your branch and open a pull request on GitHub.
