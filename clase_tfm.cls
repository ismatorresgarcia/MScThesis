%%% Versión V1.2 de la clase .cls
%%% Plantilla para el TFM de Ismael Torres García en la UPM-ETSII
\NeedsTeXFormat{LaTeX2e}
\ProvidesClass{clase_tfm}[2023/04/21 Clase TFM, V1.2, Ismael Torres García]

%%% clase de base (KOMA-Script)
\LoadClass{scrbook}                   %% clase libro KOMA-Script para el documento base
\KOMAoptions{     
    paper=a4,                         %% tamaño del papel
    fontsize=12pt,                    %% tamaño de la fuente
    twoside=true,                     %% recto-verso
    parskip=half,                     %% sangría y espaciado de parrafos
    captions=tableheading,            %% descripción de tablas arriba
    toc=listof,                       %% figuras y tablas en el ToC
    toc=bibliographynumbered,         %% bibliografía como capítulo numerado en el ToC
    chapterprefix=false,              %% nombres de los capítulos (por defecto en scrreport)
    %numbers=noenddot                 %% números de secciones sin punto
}
%%% paquetes basicos
%% fórmulas, símbolos y teoremas matemáticos.'unicode-math' sobreescribirá las definiciones necesarias después
\RequirePackage{amsmath}              %% fórmulas y entornos matemáticos. Cargar antes de 'unicode-math' por recomendación de este  
\RequirePackage{amssymb}              %% fuentes y símbolos matemáticos

\RequirePackage{unicode-math}         %% gestión de matemáticas. Los paquetes de la AMS son incompatibles, así que declaralos antes para que sobreescriba lo que considere necesario. También carga 'fontspec' para gestionar las fuentes de texto con XeLaTeX o LuaLaTeX
\unimathsetup{                        %% las opciones del paquete 'unicode-math' siguiendo la norma ISO
    math-style=ISO, 
    bold-style=ISO
} 
\setmainfont{Minion Pro}              %% Fuente comercial MinionPro de Adobe para el texto
% fuente matemática del documento
\RequirePackage{unicode-minionmath}   %% Paquete comercial MinionMath de Johannes Küster con comandos adicionales en XeLaTeX/LuaLaTeX
\setmathfont{MinionMath-Regular}[     %% Fuente comercial MinionMath de Johannes Küsner para las matemáticas
    Extension = .otf,
    Path = ./MinionMath/,
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
\setmathfont{MinionMath-Bold} [                         
    Extension = .otf,
    Path = ./MinionMath/,
    version = bold
]
\setmathfont{STIX Two Math}[          %% Fuente STIX Two Math para los estilos matemáticos que faltan en MinionMath
    range = {\symcal,\symfrak,\symscr,\symbfcal,\symbffrak,\symbfscr}
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
\RequirePackage{scrhack}              %% macros personalizadas (KOMA-Script)
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
    top=25mm,                         %% margen superior
    bottom=25mm,                      %% margen inferior
    inner=25mm,                       %% margen interior
    outer=25mm,                       %% margen exterior
    headsep=5mm,                      %% separación entre texto y encabezado
    footskip=10mm                     %% separación entre texto y pie de pagina
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
%% colores del documento
\definecolor{miazul}{RGB}{0,26,153}
\definecolor{miazull}{RGB}{119,141,255}
\definecolor{mirojo}{RGB}{153,0,0}
\definecolor{mirojoo}{RGB}{199,102,80}
\definecolor{mimorado}{RGB}{102,0,205}
\definecolor{mimoradoo}{RGB}{204,153,255}
\definecolor{miverde}{RGB}{0,153,0}
\definecolor{minaranja}{RGB}{252,147,3}
%% colores de códigos
% \definecolor{codverde}{rgb}{0,0.6,0}
% \definecolor{codazul}{rgb}{0,0,0.9}
% \definecolor{codgris}{rgb}{0.5,0.5,0.5}
% \definecolor{codmorado}{rgb}{0.58,0,0.82}
\definecolor{colorfondo}{rgb}{0.95,0.95,0.92}

%%% creación de macros y paquetes personalizados como entornos, por ejemplo
\RequirePackage[export]{adjustbox}          %% para definir un tamaño mínimo o máximo
\RequirePackage{xpatch}                     %% recomendado por el paquete de biblatex. También extiende algunas funciones del paquete 'etoolbox'

%%% miscelánea
\RequirePackage{hologo}                     %% macro '\hologo! para escribir logos especiales (LATEX, XELATEX, LUALATEX, etc)

%%% configuración de la bibliografía con 'biblatex'. Biblatex carga 'etoolbox' entre otros paquetes
\RequirePackage[                            %% gestion de la bibliografía con biblatex
    backend=biber,                          %% compilador por defecto con biblatex (recomendado)
    natbib=true,                            %% comandos del paquete natbib
    maxnames=10,                            %% número max. de nombres de autores
    minnames=1,                             %% número mix. de nombres de autores
    sorting=nyt,                            %% ordenar por nombre, año, título
    autocite=superscript,                   %% estilos de las citas del texto como superíndices
    citestyle=numeric-comp,                 %% estilo de las citas numérico tipo [1-3] (pon numeric-comp, en draft mientras edito)
    bibstyle=numeric,                       %% estilo de bibliografía numérica (pon numeric, en draft mientras edito)  
    backref=false                           %% sacar las páginas de las citas en la bibliografía (cf. p.)
]{biblatex}                   
\addbibresource{Referencias/biblio.bib}     %% archivo bibliografico

%%% matemáticas, física y química
\RequirePackage{siunitx}                    %% unidades del sistema internacional
\sisetup{                                   %% opciones de las unidades
    output-decimal-marker = {,},
    per-mode = symbol
}
\RequirePackage[version=4]{mhchem}          %% escribir reacciones, compuestos químicos, etc.

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
      frame=lines,
      framesep=2mm,
      baselinestretch=1.0,
      bgcolor=colorfondo,
      fontsize=\footnotesize,
      linenos,
}
%%% cambio el nombre que usa 'caption' para la leyenda y el índice con la opción 'newfloat' de 'minted'
\SetupFloatingEnvironment{listing}{name=Código}    
\SetupFloatingEnvironment{listing}{listname=Índice de códigos}
%% un entorno para sustituir a 'listings' cuando el código ocupa varias páginas
% \newenvironment*{code}{\captionsetup{type=listing}}{}
  
%%% PDF
\RequirePackage{pdfpages}                   %% insertar archivos pdf
\RequirePackage{hyperref}                   %% pdf con hipervínculos, referencias cruzadas, etc.
\hypersetup{
    %backref=true,                          %% añadir hipervínculos (por defecto)
    %pagebackref=true,                      %% en la bibliografía (por defecto)
    %hyperindex=true,                       %% en el índice (por defecto)
    %bookmarks=true,                        %% marcadores de Acrobat (por defecto)
    breaklinks=true,                        %% romper línea si es demasiado largo el hipervínculo
    colorlinks=true,                        %% color de hipervínculos
    urlcolor=minaranja,                     %% color de hipervínculos de direcciones web
    citecolor=miverde,                      %% color de citas bibliográficas internas
    linkcolor=mimorado,                     %% color de vínculos internos
    filecolor=mirojo,                       %% color de urls que abren archivos
    menucolor=mirojoo,                      %% color de las pestañas en Adobe
    anchorcolor=black,                      %% color del texto anclado cuando pasas por encima
    bookmarksopen=false,                    %% las pestañas de los marcadores del pdf aparecen todas desplegadas
    linktocpage=false,                      %% hipervínculos en el número de página en lugar del texto del TOC
%% metacampos del pdf :                     %% ATENCIÓN: hay que completarlos 
    pdftitle    = {Modelización de armónicos de alto orden en plasmas},
    pdfauthor   = {Ismael Torres García},
    pdfsubject  = {Trabajo Fin de Master},
    pdfkeywords = {SXRL, Kriptón, XUV, HOH, HHG, Dagon},
}
\RequirePackage{bookmark}                   %% gestión moderna de los marcadores del paquete 'hyperref'

%%% configuración de los idiomas con 'polyglossia'
\setdefaultlanguage[
    variant=spanish,
    spanishoperators=all,
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
\NewDocumentCommand\appendixprefixintoc{}             %% para añadir la palabra 'Anexo' antes de las letras del apéndice en el ToC
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
\addtokomafont{disposition}{\normalfont\bfseries}

\addtokomafont{chapter}{\fontsize{25pt}{30pt}\selectfont\color{mimorado}}
\newkomafont{chapternumber}{\fontsize{60}{72}\selectfont\color{mimoradoo}}
\setkomafont{descriptionlabel}{\normalfont\bfseries} %% Cambio el estilo de los argumentos en itemize y description, 
%% que es lo que utiliza el paquete 'glossaries' también. Por defecto Koma-Script utiliza \sffamily\bfseries
% \addtokomafont{chapterentry}{\sffamily}
% \newkomafont{chaptername}{\itshape\sffamily\normalsize\color{migris}}
    
\addtokomafont{section}{\fontsize{17pt}{21pt}\selectfont\color{mimorado}}
% \newkomafont{sectionnumber}{\fontsize{18pt}{22pt}\selectfont\color{white}}

\addtokomafont{subsection}{\fontsize{14pt}{16pt}\selectfont\color{mimorado}}
% \newkomafont{subsectionnumber}{\fontsize{14pt}{16pt}\selectfont\color{white}}

\addtokomafont{subsubsection}{\fontsize{11pt}{13pt}\selectfont\color{mimorado}}
% \newkomafont{subsubsectionnumber}{\fontsize{14pt}{16pt}\selectfont\color{white}}

\addtokomafont{paragraph}{\selectfont\color{mimorado}}
\addtokomafont{subparagraph}{\selectfont\color{mimorado}}

%%% comando para configurar el número o letra de cada capítulo en los títulos con KOMA-Script
\renewcommand*{\chapterformat}{%
  {
    \colorbox{white}{
     \parbox[c][60pt]{70pt}{\centering
      {\usekomafont{chapternumber}{\thechapter}
        }%
      }%
    }%
  }%
  \enskip
}%

%%% comando para configurar el nombre de los títulos de los capítulos con KOMA-Script
\newbox{\chapternumberbox}
\makeatletter
\renewcommand*{\chapterlinesformat}[3]{%
  \sbox\chapternumberbox{#2}{%
    \@hangfrom{\usebox\chapternumberbox}{%
      \parbox[c]{\dimexpr\textwidth-\wd\chapternumberbox\relax}{\raggedchapter\MakeUppercase{#3}}
    }%
  }%
}
\makeatother

%%% comandos para modificar el estilo de los encabezados y pies de página con KOMA-Script
\setkomafont{pageheadfoot}{\normalfont\small\bfseries\color{mimorado}}
\setkomafont{pagenumber}{\normalfont\small\bfseries\color{mimorado}}
\addtokomafont{headsepline}{\color{black}}
\addtokomafont{footsepline}{\color{black}}

%%% comandos para las páginas del título y el resumen con KOMA-Script
\newcommand*{\HRule}{\rule{\linewidth}{1pt}}
\setkomafont{subject}{\Large}
\setkomafont{dedication}{\Large\itshape}

%%% comandos para modificar nombres de macros en español con 'polyglossia'
\gappto\captionsspanish{%
  \renewcommand*{\appendixname}{Anexo}%
  \renewcommand*{\tablename}{Tabla}%
  \renewcommand*{\listtablename}{Índice de tablas}%
  % \renewcommand*{\lstlistingname}{Código}% 
  % \renewcommand*{\lstlistlistingname}{Índice de códigos}% 
}
\renewcommand*{\figureformat}{\figurename~\thefigure\autodot}
\renewcommand*{\tableformat}{\tablename~\thetable\autodot}
%\renewcommand*{\captionformat}{~--~}
\addtokomafont{caption}{\small}
\setkomafont{captionlabel}{\bfseries}

%%% comandos para el espacio entre párrafos y secciones con KOMA-Script
\RedeclareSectionCommand[%
  tocnumwidth=1.75em,
  beforeskip=-.25\baselineskip ,
  afterskip=2.5\baselineskip plus .1\baselineskip minus 
  .12\baselineskip
]{chapter}

\RedeclareSectionCommands[
  tocindent=1.75em,
  tocnumwidth=2.5em,
  beforeskip=-.75\parskip plus -1ex minus -.2ex,
  afterskip=.5\parskip plus .2ex
]{section}

\RedeclareSectionCommands[
  tocindent=4.25em,
  tocnumwidth=3.5em,
  beforeskip=-.75\parskip plus -1ex minus -.2ex,
  afterskip=.5\parskip plus .2ex
]{subsection}

\RedeclareSectionCommands[
  beforeskip=1sp plus -1sp minus 1sp,
  afterskip=1sp plus -1sp minus 1sp,
]{paragraph}

%%% Definición de comandos útiles para escribir
%% conjuntos de números
\newcommand*{\C}{\symbb{C}}
\newcommand*{\K}{\symbb{K}}
\newcommand*{\N}{\symbb{N}}
\newcommand*{\Q}{\symbb{Q}}
\newcommand*{\R}{\symbb{R}}
\newcommand*{\Z}{\symbb{Z}}
\newcommand*{\X}{\symbb{X}}
%% números complejos
\DeclareMathOperator{\RE}{\operatorname{Re}}
\DeclareMathOperator{\IM}{\operatorname{Im}}
\newcommand*{\iu}{\symrm{i}}
\newcommand*{\eu}{\symrm{e}}
%% operadores matemáticos
\newcommand*{\dif}{\symrm{d}}
\newcommand*{\Dif}{\symrm{D}}
\newcommand*{\del}{\symrm{\Delta}}
\newcommand*{\nab}{\symbf{\nabla}}




