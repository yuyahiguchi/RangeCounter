#!/bin/sh

if [ $# -lt 3 ];then
    echo "usage1 : ./RCscript.sh [absorberlength] [degraderlength] [material] [pid]"
    echo "This is used to run for a paticular particle.\nNote that pid means PDGid (m- = 13, mu+ = -13, pi+ = 211, pi- = -211)\n"
    echo "usage2 : ./RCscript.sh [absorberlength] [degraderlength] [material]"
    echo "This is used to run mu-, mu+, pi+ and pi-.\n"
    exit 0
fi

if [ $1 == "-h" || $1 == "--help" ];then
    echo "usage1 : ./RCscript.sh [absorberlength] [degraderlength] [material] [pid]"
    echo "This is used to run for a paticular particle\n\n"
    echo "usage2 : ./RCscript.sh [absorberlength] [degraderlength] [material]"
    echo "This is used to run mu-, mu+, pi+ and pi-.\n"
    exit 0
fi

if [ $4 == "13" ];then
    particle="mu-"
elif [ $4 == "-13" ];then
    particle="mu+"
elif [ $4 == "211" ];then
    particle="pi+"
elif [ $4 == "-211" ];then
    particle="pi-"
fi


SimPath="/Users/higuchiyuya/range-counter/higuchi/data/Simulation"
DataPath="${SimPath}/Absorber$1"

if [ ! -d ${DataPath} ];then
    mkdir ${DataPath}
fi

if [ $# -eq 3 ];then
    g4bl Sim_for_RealMonteCarlo.in hour=1 absorberlength=$1 degraderlength=$2 material=$3 datadir=${DataPath} pid=211 
    #macro/CompressRootFile ${DataPath}/Range-counter-hour1-degrader$2\[mm\]-pi+-in$3.root
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++"
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++"
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++"
    
    g4bl Sim_for_RealMonteCarlo.in hour=1 absorberlength=$1 degraderlength=$2 material=$3 datadir=${DataPath} pid=-13 
    #macro/CompressRootFile ${DataPath}/Range-counter-hour1-degrader$2\[mm\]-mu+-in$3.root
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++"
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++"
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++"

    g4bl Sim_for_RealMonteCarlo.in hour=1 absorberlength=$1 degraderlength=$2 material=$3 datadir=${DataPath} pid=13
    #macro/CompressRootFile ${DataPath}/Range-counter-hour1-degrader$2\[mm\]-mu--in$3.root
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++"
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++"
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++"
    
    g4bl Sim_for_RealMonteCarlo.in hour=1 absorberlength=$1 degraderlength=$2 material=$3 datadir=${DataPath} pid=-211 
    #macro/CompressRootFile ${DataPath}/Range-counter-hour1-degrader$2\[mm\]-pi--in$3.root
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++"
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++"
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++"

elif [ $# -eq 4 ];then
    g4bl Sim_for_RealMonteCarlo.in hour=1 absorberlength=$1 degraderlength=$2 material=$3 datadir=${DataPath} pid=$4
    #macro/CompressRootFile ${DataPath}/Range-counter-hour1-degrader$2\[mm\]-${particle}-in$3.root
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++"
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++"
    echo "+++++++++++++++++++++++++++++++++++++++++++++++++"

fi
