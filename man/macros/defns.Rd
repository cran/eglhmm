\newcommand{\today}{\Sexpr[stage=build]{format(Sys.Date(),"\\\%d/\\\%m/\\\%Y")}}
\newcommand{\icsl}{\out{\hspace*{1em}}}
\newcommand{\icsh}{\out{<a>&emsp; &emsp;</a>}}
\newcommand{\ics}{\ifelse{latex}{\icsl}{\ifelse{html}{\icsh}{ }}}
