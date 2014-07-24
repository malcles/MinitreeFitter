#!/bin/bash


usage() {
    echo "`basename $0` --help"
}

inDir=/afs/cern.ch/work/f/fcouderc/public/diphoton2012_cicpfSel_cms53x_v1.0/
outDir=FitIncl/Spin2pm_mu

interactive=0

if ! options=$( getopt -o o:i -l help -- "$@" )
then
    # something went wrong, getopt will put out an error message for us
    exit 1
fi

eval set -- $options
while [ $# -gt 0 ]; do
    case "$1" in
	-h | --help) usage; exit 0;;
	-o) outDir=$2; shift;;
	-i) interactive=1;;
	(--) shift; break;;
	(-*) echo "$0: error - unrecognized option $1" 1>&2; usage >> /dev/stderr; exit 1;;
	(*)  break;;
    esac
    shift
done

###### print and check configuration
echo "================================="
echo " output dir: "$outDir


optionsub="-q 1nd "

dirafs=`pwd`
dirscript=tmpBatchOut/
eosPref=root://eoscms//eos/cms
castorPref=


config() {
    mu=$1
    script=$2
    outD=${outDir}${mu}
    inD=${inDir}
    options="--inDir=${inDir} --simFit=1 --sigHyp=1 --scaleSig=${mu}  --outDir=${outDir}${mu} -i 0"
    exe="MiniTreeFitterSpin ${options}"
    echo "$exe"

    cat > $script<<EOF 
#!/bin/bash
cd $dirafs
source etc/scripts/setup.sh
gcccom="\`which gcc\`"
echo "gcc:" \$gcccom
echo "where am I:\`pwd\`"
echo  $exe
$exe
EOF
    chmod +x $script
}


mkdir -p $dirscript

for mu in $( seq 1.0 0.1 1.0 ); do 
    script=script$$
    script=${script}_${spin2InclFit}_${mu}
    cd $dirscript/
    config $mu $script
    echo "-> created script: "$script
    if [ $interactive -eq 1 ]; then
	source $script
    else
	bsub $optionsub $script 
    fi
    cd -
done


