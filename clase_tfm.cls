%%% Versión V1.2 de la clase .cls
%%% Plantilla para el TFM de Ismael Torres García en la UPM-ETSII
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{clase_tfm}[2023/04/21 Clase TFM, V1.2, Ismael Torres García]

%%% clase de base (KOMA-Script)
\LoadClass{scrbook}                   %% clase libro KOMA-Script para el documento base
\KOMAoptions{     
  paper = a4,                         %% tamaño del papel
  fontsize = 12pt,                    %% tamaño de la fuente
  twoside = true,                     %% recto-verso
  parskip = half,                     %% sangría y espaciado de párrafos
  chapterprefix = false,              %% nombres de los capítulos (por defecto en scrreport)
  captions = tableheading,            %% descripción de tablas arriba
  toc = listof,                       %% figuras y tablas en el ToC
  toc = bibliographynumbered,         %% bibliografía como capítulo numerado en el ToC
  numbers = noenddot                  %% números de capítulos, secciones, ... sin punto
}
%%% paquetes básicos
%% fórmulas, símbolos y teoremas matemáticos.'unicode-math' sobrescribirá las definiciones necesarias después
\RequirePackage{amsmath}              %% fórmulas y entornos matemáticos. Cargar antes de 'unicode-math' por recomendación de este  
\RequirePackage{amssymb}              %% fuentes y símbolos matemáticos (carga 'amsfonts' también)

\RequirePackage{unicode-math}         %% gestión de matemáticas. Los paquetes de la AMS son incompatibles, así que ponlos antes para que sobrescriba lo que considere necesario. También carga 'fontspec' para gestionar las fuentes de texto con XeLaTeX o LuaLaTeX
\unimathsetup{                        %% las opciones del paquete 'unicode-math' siguiendo la norma ISO
  math-style = ISO, 
  bold-style = ISO
} 
\setmainfont{Minion Pro}[             %% Fuente comercial Minion Pro de Adobe (desarrollador Robert Slimbach) para el texto en Serif
  Renderer = OpenType
]
\setsansfont{Cronos Pro}[             %% Fuente comercial Cronos Pro de Adobe (desarrollador Robert Slimbach) para el texto en Sans Serif
  Renderer = OpenType
]              
\RequirePackage{unicode-minionmath}   %% Paquete comercial Minion Math de Typoma (desarrollador Johannes Küster) con comandos adicionales en XeLaTeX/LuaLaTeX
%% fuente matemática del documento
\setmathfont{MinionMath-Regular}[     %% Fuente comercial Minion Math de Johannes Küsner para las matemáticas
  Extension = .otf,
  Path = /Users/ytoga/Library/Fonts/,
  Scale = 1,
  Script = Math,
  SizeFeatures = {
    {Size = -6, Font = MinionMath-Tiny,
    Style = MathScriptScript},
    {Size = 6-8.4, Font = MinionMath-Capt,
    Style = MathScript},
    {Size = 8.4-, Font = MinionMath-Regular}
  }
]
\setmathfont{MinionMath-Bold}[                         
  Extension = .otf,
  Path = /Users/ytoga/Library/Fonts/,
  version = bold
]
\setmathfont{STIXTwoMath}[          %% Fuente STIX Two Math para los estilos matemáticos que faltan en MinionMath
  Extension = .otf,
  Path = /System/Library/Fonts/Supplemental/,
  range = {\symcal,\symfrak,\symbfcal,\symbffrak},
  StylisticSet = 1
]

%%% gestión de los idiomas. Las opciones de 'polyglossia' las cargo más adelante por recomendación del manual
\RequirePackage{polyglossia}

%%% microtipografía de las fuentes con 'microtype'. Con LuaLaTeX tiene más funcionalidades (expansión y protrusión) que XeLaTeX (solo protrusión)
\RequirePackage{microtype}            %% microtipografía (detalles finos entre palabras, en los márgenes, etc.)

%%% comandos para introducir datos importantes
\newcommand*{\nomautor}{Ismael Torres García}
\newcommand*{\nomtutor}{Eduardo Oliva Gonzalo}
\newcommand*{\dirinst}{C/ José Gutiérrez Abascal, 2}
\newcommand*{\cpinst}{28006 Madrid (España)}
\newcommand*{\titproy}{Trabajo Fin de Master}
\newcommand*{\titenc}{Estudio del impacto de la forma de un canal de plasma en la amplificación de un haz XUV}
\newcommand*{\unienc}{Universidad Politécnica de Madrid}
\newcommand*{\etsiienc}{Escuela Técnica Superior de Ingenieros Industriales (UPM)}
\newcommand*{\AUA}{Abreviaturas, unidades y acrónimos}

%%% configuración y disposición del texto con KOMA-Script
\RequirePackage{scrhack}              %% macros especiales y personalizadas (KOMA-Script)
\RequirePackage{scrdate}              %% comandos para las fechas (KOMA-Script)
\RequirePackage[  
  % automark,                       %% encabezados automáticos \automark[section]{chapter} para scrbook
  % autooneside=false,              %% desactivo esta opción para que automark controle los encabezados de capítulos y secciones
  headsepline,                      %% línea horizontal del encabezado
  footsepline                       %% línea horizontal del pie de página
]{scrlayer-scrpage}                   %% encabezados y pies de página personalizados (KOMA-Script)
\pagestyle{scrheadings}               %% estilo personalizado para todas las páginas menos la primera de un capítulo (que utiliza plain.scrheadings)
\automark[chapter]{chapter}           %% nombre de capítulos en encabezados pares e impares
\automark*[]{section}                 %% nombre de secciones en encabezados cuando aparezcan (cambio el 'leftmark' de capítulo por sección)
\lefoot[]{\pagemark}
\refoot[]{\etsiienc}
\rohead[]{\titenc}
\lofoot[]{\nomautor}
\rofoot[\pagemark]{\pagemark}

\RequirePackage[
  top = 25mm,                         %% margen superior
  bottom = 25mm,                      %% margen inferior
  inner = 25mm,                       %% margen interior
  outer = 25mm,                       %% margen exterior
  headsep = 5mm,                      %% separación entre texto y encabezado
  footskip = 10mm                     %% separación entre texto y pie de pagina
]{geometry}                           %% geometría de las páginas

%%% tablas y figuras
\RequirePackage{graphicx}             %% imágenes y figuras
\graphicspath{{Figuras/}}             %% ruta de las figuras
\RequirePackage{caption}              %% leyenda de las figuras. Está cargado para utilizar 'subcaption' después con la última versión
\RequirePackage{subcaption}           %% varias subfiguras
\RequirePackage{booktabs}             %% tablas formales
\RequirePackage{multirow}             %% fusión de varias filas en una celda
\RequirePackage{longtable}            %% tablas muy largas
\RequirePackage{array}                %% tablas
\RequirePackage{booktabs}             %% tablas profesionales
\RequirePackage[table]{xcolor}        %% mejores colores. 'table' para colorear filas y columnas de tablas
\RequirePackage{tikz}                 %% para hacer dibujos (aquí solo lo utilizo para los capítulos)
%% colores del documento
\definecolor{miazul}{RGB}{0,30,155}
\definecolor{miazulll}{RGB}{0,50,85}
\definecolor{mirojo}{RGB}{153,0,0}
\definecolor{mirojoo}{RGB}{200,100,80}
\definecolor{miverde}{RGB}{0,153,0}
\definecolor{mimorado}{RGB}{140,45,165}
\definecolor{minaranja}{RGB}{250,150,5}
\colorlet{miazull}{white!85!miazul}
\colorlet{secboxcolor}{miazul}
\colorlet{secnumcolor}{white}
%% colores de códigos
% \definecolor{codverde}{rgb}{0,0.6,0}
% \definecolor{codazul}{rgb}{0,0,0.9}
% \definecolor{codgris}{rgb}{0.5,0.5,0.5}
% \definecolor{codmorado}{rgb}{0.58,0,0.82}
\definecolor{colorfondo}{rgb}{0.95,0.95,0.92}

%%% creación de macros y paquetes personalizados como entornos, por ejemplo
\RequirePackage[export]{adjustbox}          %% para definir un tamaño mínimo o máximo en las cajas
\RequirePackage{xpatch}                     %% recomendado por el paquete de 'biblatex'. También extiende algunas funciones del paquete 'etoolbox'

%%% más personalización
\RequirePackage[lines=3]{lettrine}                   %% para escribir los comienzos de capítulos en letras grandes
\renewcommand*{\LettrineFontHook}{\fontspec{Cronos Pro}\bfseries\color{miazul}} %% para el estilo de la letra inicial
\renewcommand*{\LettrineTextFont}{\sffamily\bfseries\fontsize{10pt}{12pt}\selectfont\color{miazul}\MakeUppercase} %% para el estilo del segundo argumento

%%% configuración de la bibliografía con 'biblatex'. Biblatex carga 'etoolbox' entre otros paquetes
\RequirePackage[                            %% gestion de la bibliografía con biblatex
  backend = biber,                          %% compilador por defecto con biblatex (recomendado)
  natbib = true,                            %% comandos del paquete natbib (autocite, por ejemplo)
  maxnames = 10,                            %% número max. de nombres de autores
  minnames = 1,                             %% número mix. de nombres de autores
  sorting = nyt,                            %% ordenar por nombre, año, título
  autocite = superscript,                   %% estilos de las citas del texto como superíndices
  citestyle = numeric-comp,                 %% estilo de las citas numérico tipo [1-3] (pon numeric-comp, en draft mientras edito)
  bibstyle = numeric,                       %% estilo de bibliografía numérica (pon numeric, en draft mientras edito)  
  backref = false                           %% sacar las páginas de las citas en la bibliografía (cf. p.)
]{biblatex}                   
\addbibresource{Zotero.bib}                 %% archivo bibligráfico que estoy utilizando

%%% matemáticas, física y química
\RequirePackage{siunitx}                    %% unidades del sistema internacional
\sisetup{                                   %% opciones de las unidades del paquete
  output-decimal-marker = {,},
  uncertainty-mode = separate,
  per-mode = symbol
}
\RequirePackage[version=4]{mhchem}          %% para escribir reacciones, compuestos químicos, etc.

%%% Algoritmos y códigos de lenguajes de programación con 'minted' o 'listings'
%% en 'minted' los entornos flotantes se invocan con \begin{listing} \end{listing} y 
%% el código dentro del entorno \begin{minted}{lenguaje} \end{minted}
%% con \inputminted{opciones}{lenguaje}{archivo.extensión} puedes introducir el código de programas guardados
%% los comandos para el nombre de leyendas y el índice de códigos son \listingscaption y \listloflistingscaption
%% (solo si la opción 'newfloat' no está cargada). 'minted' utiliza una librería externa llamada 'Pygment' para
%% resaltar la sintaxis (un programa de python)
\RequirePackage[newfloat]{minted}             %% 'newfloat' es parte del paquete 'caption'
\usemintedstyle{emacs}                        %% estilo predeterminado de 'minted' con 'Pygment' 
\setminted{                                   %% opciones globales para el formato del código
  frame = lines,
  framesep = 2mm,
  baselinestretch = 1.0,
  bgcolor = colorfondo,
  fontsize = \footnotesize,
  linenos,
}
%%% cambio el nombre que usa 'caption' para la leyenda y el índice con la opción 'newfloat' de 'minted'
\SetupFloatingEnvironment{listing}{name=Código}    
\SetupFloatingEnvironment{listing}{listname=Índice de códigos}
  
%%% PDF
\RequirePackage{pdfpages}                   %% insertar archivos pdf
\RequirePackage{bookmark}                   %% gestión moderna de los marcadores del paquete 'hyperref'
\RequirePackage{hyperref}                   %% pdf con hipervínculos, referencias cruzadas, etc.
\hypersetup{
  %backref=true,                          %% añadir hipervínculos (por defecto)
  %pagebackref=true,                      %% en la bibliografía (por defecto)
  %hyperindex=true,                       %% en el índice (por defecto)
  %bookmarks=true,                        %% marcadores de Acrobat (por defecto)
  breaklinks = true,                        %% romper línea si es demasiado largo el hipervínculo
  colorlinks = true,                        %% color de hipervínculos
  urlcolor = minaranja,                     %% color de hipervínculos de direcciones web
  citecolor = miverde,                      %% color de citas bibliográficas internas
  linkcolor = miazul,                       %% color de vínculos internos
  filecolor = mirojo,                       %% color de urls que abren archivos
  menucolor = mirojoo,                      %% color de las pestañas en Adobe
  anchorcolor = black,                      %% color del texto anclado cuando pasas por encima
  bookmarksopen = false,                    %% las pestañas de los marcadores del pdf aparecen todas desplegadas
  %linktocpage=false,                     %% hipervínculos en el número de página en lugar del texto del TOC
%% metacampos del pdf :                     %% ATENCIÓN: hay que completarlos 
  pdftitle    = {Modelización de armónicos de alto orden en plasmas},
  pdfauthor   = {Ismael Torres García},
  pdfsubject  = {Trabajo Fin de Master},
  pdfkeywords = {SXRL, Kriptón, XUV, HOH, HHG, Dagon},
}
\AfterTOCHead[toc]{\hypersetup{linkcolor=.}}%% ejecuta el comando \hypersetup justo después del ToC (para activar los links, pero utilizar el color azul de linkcolor solo en los capítulos) con KOMA-Script

%%% configuración de los idiomas con 'polyglossia'
\setdefaultlanguage[
  variant = spanish,
  spanishoperators = all,
]{spanish}
\setotherlanguage[variant=british]{english}
\RequirePackage[autostyle]{csquotes}        %% entorno para comillas tipográficas (y más cosas) 

%%% Lista de acrónimos, constantes y unidades con 'glossaries'
\RequirePackage[
  toc,
  acronym,
  nonumberlist,
  subentrycounter
]{glossaries}
\makeglossaries
\setacronymstyle{long-short}
\loadglsentries{Preliminares/acronimos.tex}

%%% comandos para configurar la entrada 'anexos' en el ToC con KOMA-Script
\newcommand*{\useprefix}[1]{#1}                       %% defino un comando para almacenar la palabra 'anexo'
\newcommand*{\appendixprefixintoc}                    %% para añadir la palabra 'Anexo' antes de las letras del apéndice en el ToC
{%
  \DeclareTOCStyleEntry[%                             %% defino el estilo de los capítulos del anexo en el ToC
    entrynumberformat=\useprefix{\appendixname~},
    dynnumwidth                                     %% mínimo ancho por defecto {default} de la entrada de la letra 'A, B, C, ...'
  ]{default}{chapter}
  \DeclareTOCStyleEntries[%                           %% defino el estilo de las secciones, subsecciones, ... en el ToC
    dynindent                                         %% mínima distancia horizontal por defecto {default} entre margen izquierdo y entradas
  ]{default}{section,subsection,paragraph}
}
\newcommand*{\appendixmore}{%                         %% introduzco la entrada 'anexo' en los nombres del ToC
  \addtocontents{toc}{\appendixprefixintoc}%
}

%%% comandos para configurar tamaños y fuentes de los distintos niveles
\addtokomafont{disposition}{\sffamily\bfseries\color{miazul}}

\renewcommand*{\raggedchapter}{\raggedleft} %% Justificación a la derecha del título de los capítulos
\setkomafont{chapter}{\sffamily\fontsize{25pt}{30pt}\selectfont\color{miazul}}
\setkomafont{chapterprefix}{\sffamily\fontsize{11pt}{13pt}\selectfont\color{miazul}}
\newkomafont{chapternumber}{\fontsize{50pt}{60pt}\selectfont\color{miazul}}

\setkomafont{descriptionlabel}{\normalfont\bfseries} %% Cambio el estilo de los argumentos en itemize y description, 
%% que es lo que utiliza el paquete 'glossaries' también. Por defecto Koma-Script utiliza \sffamily\bfseries
    
\addtokomafont{section}{\sffamily\fontsize{17pt}{21pt}\selectfont\color{miazul}}
\addtokomafont{subsection}{\sffamily\fontsize{14pt}{17pt}\selectfont\color{miazul}}
\addtokomafont{subsubsection}{\sffamily\fontsize{11pt}{13pt}\selectfont\color{miazul}}
\addtokomafont{paragraph}{\normalfont\sffamily\itshape\fontsize{12pt}{14pt}\selectfont\color{miazul}}
\addtokomafont{subparagraph}{\normalfont\sffamily\itshape\fontsize{12pt}{14pt}\selectfont\color{miazul}}

%%% empiezo a configurar con TikZ el formato de cada capítulo en los títulos 
\tikzset{
  headings/base/.style = {
    outer sep = -5pt,
    inner sep = 10pt,
  },
  headings/chapterbackground/.style = {
    headings/base,
    shade,
    left color = white,
    right color = miazull,
  },
  headings/chapapp/.style = {
    headings/base,
    font = \usekomafont{chapterprefix},
  },
  headings/chapternumber/.style = {
    headings/base,
    font = \usekomafont{chapternumber},
  },
  headings/chapterline/.style = {
    miazull,
    line width = 2.5pt,
  }
}
%%% comandos para configurar el formato de cada capítulo en los títulos con KOMA-Script
\makeatletter
\renewcommand*\chapterlinesformat[3]{%
  \Ifstr{#1}{chapter}{%
    \begin{tikzpicture}[baseline=(title.base)]
      \node[headings/chapterbackground](title){%
        \parbox[c][\height]
          {\dimexpr\textwidth-1\pgfkeysvalueof{/pgf/inner xsep}\relax}
          {\raggedchapter #3}%
      };
      \node[headings/chapapp,anchor=south east]
      at (title.north east){\Ifstr{#2}{}{}{\MakeUppercase{\chapapp}}\strut};
      \useasboundingbox
        (current bounding box.south west)
        rectangle
        ([yshift=-10pt]current bounding box.south east);
      \draw[headings/chapterline]
        (current bounding box.south east)++(+.5\pgflinewidth,0)--+(0,\paperheight);
      \node[anchor=base west,headings/chapternumber]
        at ([xshift=10pt]title.base-|current bounding box.east){#2};
    \end{tikzpicture}
    \par
  }{%
    \@hangfrom{#2}{#3}% los niveles de otras secciones que utilicen style=chapter
  }%
}
\makeatother

%%% comando para configurar el estilo de las secciones y subsecciones con KOMA-Script
\renewcommand*{\sectionformat}{\colorbox{secboxcolor}{\textcolor{secnumcolor}{~\thesection\autodot}~}\quad}
\renewcommand*{\subsectionformat}{\colorbox{secboxcolor}{\textcolor{secnumcolor}{~\thesubsection\autodot}~}\quad}

\newsavebox{\secnumbox}
\newlength{\secrulewidth}
\setlength{\secrulewidth}{1pt}
\newcommand*{\boxedandruledsec}[2]{%
  \IfUseNumber{\leavevmode\rlap{\color{secboxcolor}\rule[\dimexpr-\fboxsep-\secrulewidth\relax]{\textwidth}{\secrulewidth}}}{}%
  \savebox{\secnumbox}{#1}%
  \parbox[b]{\wd\secnumbox}{\usebox{\secnumbox}}%
  \parbox[b]{\dimexpr\textwidth-\wd\secnumbox\relax}{#2}\par\nobreak%
}

\newcommand\originalsectionlinesformat{}
\let\originalsectionlinesformat\sectionlinesformat

\renewcommand*{\sectionlinesformat}[4]{\setlength\fboxsep{\dp\strutbox}%
  \Ifstr{#1}{section}{\boxedandruledsec{#3}{#4}}{%
    \Ifstr{#1}{subsection}{\boxedandruledsec{#3}{#4}}{%
      \originalsectionlinesformat{#1}{#2}{#3}{#4}% los niveles de otras secciones
    }}%
}

%%% comandos para modificar el estilo de los encabezados y pies de página con KOMA-Script
\setkomafont{pageheadfoot}{\sffamily\small\bfseries\selectfont\color{miazul}}
\setkomafont{pagenumber}{\sffamily\small\bfseries\selectfont\color{miazul}}
\addtokomafont{headsepline}{\selectfont\color{miazul}}
\addtokomafont{footsepline}{\selectfont\color{miazul}}

%%% comandos para las páginas del título y el resumen con KOMA-Script
\setkomafont{subject}{\Large}
\setkomafont{dedication}{\Large\itshape}
\setkomafont{title}{\normalfont\bfseries} %% Cambio el estilo del title en la portada
\setkomafont{titlehead}{\normalsize} %% Cambio el estilo del titlehead en la portada

%%% comandos para modificar los nombres de descripciones en español con 'polyglossia'
\gappto\captionsspanish{%
  \renewcommand*{\appendixname}{Anexo}%
  \renewcommand*{\tablename}{Tabla}%
  \renewcommand*{\listtablename}{Índice de tablas}%
}

\renewcommand*{\figureformat}{\figurename~\thefigure\autodot}
\renewcommand*{\tableformat}{\tablename~\thetable\autodot}
\renewcommand*{\captionformat}{~\quad~}
\addtokomafont{caption}{\footnotesize\sffamily\color{miazulll}}
\setkomafont{captionlabel}{\footnotesize\sffamily\bfseries\color{mimorado}}
\setcapindent{0pt} %% Pongo a 0 el sangrado del cuerpo de las figuras y tablas con Koma-Script

%%% comandos para el espacio entre párrafos y secciones con KOMA-Script
\RedeclareSectionCommand[%
  tocnumwidth = 1.75em,
  beforeskip = -.25\baselineskip,
  afterskip = 2.5\baselineskip plus .1\baselineskip minus .12\baselineskip
]{chapter}

\RedeclareSectionCommands[
  tocentryformat = \small\sffamily,
  toclinefill = \hfill, % me cargo la línea de puntos del ToC
  tocindent = 1.75em,
  tocnumwidth = 2.5em,
  beforeskip = -1.0\parskip plus -1ex minus -.2ex,
  afterskip = 1sp plus -1sp minus 1sp
]{section}

\RedeclareSectionCommands[
  tocentryformat = \small\sffamily,
  toclinefill = \hfill, % me cargo la línea de puntos del ToC
  tocindent = 4.25em,
  tocnumwidth = 3.5em,
  beforeskip = -1.0\parskip plus -1ex minus -.2ex,
  afterskip = 1sp plus -1sp minus 1sp
]{subsection}

\RedeclareSectionCommands[
  beforeskip = -1.0\parskip plus -1ex minus -.2ex,
  afterskip = 1sp plus -1sp minus 1sp
]{paragraph}

%%% Definición de comandos útiles para escribir con los snippets
%% vectores de la base
\newcommand*{\unitvector}{\symbf{e}}
\newcommand*{\uvec}{\hat{e}}
%% derivadas totales
\newcommand*{\dvOne}[1]{\frac{\operatorname{d}}{\operatorname{d}\!#1}}
\newcommand*{\dv}[2]{\frac{\operatorname{d}\!#1}{\operatorname{d}\!#2}}
\newcommand*{\dvN}[3]{\frac{\operatorname{d^{#3}}\!#1}{\operatorname{d}\!{#2}^{#3}}}
%% derivadas parciales
\newcommand*{\pdvOne}[1]{\frac{\partial}{\partial#1}#1}
\newcommand*{\pdv}[2]{\frac{\partial#1}{\partial#2}}
\newcommand*{\pdvN}[3]{\frac{\partial^{#3}#1}{\partial{#2}^{#3}}}
%% diferenciales de primer orden
\newcommand*{\diff}{\operatorname{d}\!}
\newcommand*{\pdiff}{\partial}
%% nombres de conceptos
\newcommand*{\RE}{\operatorname{Re}}
\newcommand*{\IM}{\operatorname{Im}}
%% conjuntos de números
\newcommand*{\iu}{\symrm{i}}
\newcommand*{\eu}{\symrm{e}}
%% operadores matemáticos
\newcommand*{\grad}{\nabla}
\newcommand*{\curl}{\nabla\times}
\newcommand*{\diver}{\nabla\cdot}
\newcommand*{\laplacian}{\increment}
%% delimitadores importantes
\newcommand*{\xval}[1]{\left\langle#1\right\rangle}
\newcommand*{\abs}[1]{\left|#1\right|}
\newcommand*{\ket}[1]{\left|#1\right\rangle}
\newcommand*{\bra}[1]{\left\langle#1\right|}
