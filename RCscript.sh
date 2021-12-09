#!/bin/sh

if [ $# -lt 3 ];then
    echo "usage1 : ./RCscript.sh [pid] [absorberlength] [degraderlength] [material]"
    echo "This is used to run for a paticular particle.\nNote that pid means PDGid (m-=13, mu+=-13, pi+=211)\n"
    echo "usage2 : ./RCscript.sh [absorberlength] [degraderlength] [material]"
    echo "This is used to run mu-, mu+ and pi+.\n"
    exit 0
fi

if [ $1 == "-h" || $1 == "--help" ];then
    echo "usage1 : ./RCscript.sh [pid] [absorberlength] [degraderlength] [material]"
    echo "This is used to run for a paticular particle\n\n"
    echo "usage2 : ./RCscript.sh [absorberlength] [degraderlength] [material]"
    echo "This is used to run mu-, mu+ and pi+\n"
    exit 0
fi

if [ $1 == "13" ];then
    particle="mu-"
elif [ $1 == "-13" ];then
    particle="mu+"
elif [ $1 == "211" ];then
    particle="pi+"
fi

pwd=`pwd`
DataPath="${pwd}/data/Simulation/Absorber$1"

if [ ! -d ${DataPath} ];then
    mkdir -p ${DataPath}
fi

if [ $# -eq 3 ];then
    g4bl Sim_for_RealMonteCarlo.in pid=211 hour=1 absorberlength=$1 degraderlength=$2 material=$3 datadir=${DataPath}
    ana/CompressRootFile ${DataPath}/Range-counter-hour1-degrader$2\[mm\]-pi+-in$3.root
    rm -f ${DataPath}/Range-counter-hour1-degrader$2\[mm\]-pi+-in$3.root
    mv ${DataPath}/Range-counter-hour1-degrader$2\[mm\]-pi+-in$3_compressed.root ${DataPath}/Range-counter-hour1-degrader$2\[mm\]-pi+-in$3.root
    echo "rm -f ${DataPath}/Range-counter-hour1-degrader$2\[mm\]-mu--in$3.root"
    echo "mv ${DataPath}/Range-counter-hour1-degrader$2\[mm\]-mu--in$3_commpressed.root ${DataPath}/Range-counter-hour1-degrader$2\[mm\]-mu--in$3.root"
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++"
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++"
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++"
    
    g4bl Sim_for_RealMonteCarlo.in pid=-13 hour=1 absorberlength=$1 degraderlength=$2 material=$3 datadir=${DataPath}
    ana/CompressRootFile ${DataPath}/Range-counter-hour1-degrader$2\[mm\]-mu+-in$3.root
    rm -f ${DataPath}/Range-counter-hour1-degrader$2\[mm\]-mu+-in$3.root
    mv ${DataPath}/Range-counter-hour1-degrader$2\[mm\]-mu+-in$3_compressed.root ${DataPath}/Range-counter-hour1-degrader$2\[mm\]-mu+-in$3.root
    echo "rm -f ${DataPath}/Range-counter-hour1-degrader$2\[mm\]-mu--in$3.root"
    echo "mv ${DataPath}/Range-counter-hour1-degrader$2\[mm\]-mu--in$3_commpressed.root ${DataPath}/Range-counter-hour1-degrader$2\[mm\]-mu--in$3.root"
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++"
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++"
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++"

    g4bl Sim_for_RealMonteCarlo.in pid=13 hour=1 absorberlength=$1 degraderlength=$2 material=$3 datadir=${DataPath}
    ana/CompressRootFile ${DataPath}/Range-counter-hour1-degrader$2\[mm\]-mu--in$3.root
    rm -f ${DataPath}/Range-counter-hour1-degrader$2\[mm\]-mu--in$3.root
    mv ${DataPath}/Range-counter-hour1-degrader$2\[mm\]-mu--in$3_compressed.root ${DataPath}/Range-counter-hour1-degrader$2\[mm\]-mu--in$3.root
    echo "rm -f ${DataPath}/Range-counter-hour1-degrader$2\[mm\]-mu--in$3.root"
    echo "mv ${DataPath}/Range-counter-hour1-degrader$2\[mm\]-mu--in$3_commpressed.root ${DataPath}/Range-counter-hour1-degrader$2\[mm\]-mu--in$3.root"
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++"
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++"
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++"

elif [ $# -eq 4 ];then
    g4bl Sim_for_RealMonteCarlo.in pid=$1 hour=1 absorberlength=$2 degraderlength=$3 material=$4 datadir=${DataPath}
    ana/CompressRootFile ${DataPath}/Range-counter-hour1-degrader$3\[mm\]-${particle}-in$4.root
    rm -f ${DataPath}/Range-counter-hour1-degrader$3\[mm\]-${particle}-in$4.root
    mv ${DataPath}/Range-counter-hour1-degrader$3\[mm\]-${particle}-in$4_compressed.root ${DataPath}/Range-counter-hour1-degrader$3\[mm\]-${particle}-in$4.root
    echo "rm -f ${DataPath}/Range-counter-hour1-degrader$3\[mm\]-${particle}-in$4.root"
    echo "mv ${DataPath}/Range-counter-hour1-degrader$3\[mm\]-${particle}-in$4_compressed.root ${DataPath}/Range-counter-hour1-degrader$3\[mm\]-${particle}-in$4.root"
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++"
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++"
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++"

fi
