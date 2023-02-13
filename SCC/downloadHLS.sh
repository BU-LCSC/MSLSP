#!/bin/bash

# Download all the HLS L30 and S30 files for a single tile ID or for all the tile IDs 
# in a given file; a year can be specified as an option.
#
# version: 1.0
#
# Junchang Ju
# Jul 26, 2017: Download for an entire test site.
# Mar 16, 2018: Modified to work for a tile ID, or a file of tile IDs.
# Oct 11, 2018: Delete any partially downloaded file on interruption, and
# 		only download new files or newer files. 
# Oct 25, 2018: Fix a bug reported by Tyler Anderson so that all years of S30 is now downloaded.

tilePara=
YEAR=			# Oct 25, 2018: change to uppercase because year is also used as a variable later.	
outdir=
if [ $# = 3 ] 
then
	# $1 is "-t"
	tilePara=$2
	OUTDIR=$3
elif [ $# = 5 ] 
then
  	# $1 is "-t"
        tilePara=$2
        # $3 is "-y"
        YEAR=$4
        OUTDIR=$5
else
	echo "Example usage: $0 -t 11SPC /data/hls" >&2
  	echo "               $0 -t 11SPC -y 2017 /data/hls " >&2
  	echo "               $0 -t mytilelist -y 2017 /data/hls" >&2
     	echo "in which " >&2
     	echo "   11SPC is a 5-character tile ID " >&2
     	echo "   /data/hls is the directory where to save the data" >&2
	echo "   2017 indicates download for year 2017 only" >&2 
     	echo "Alternatively, a text file (named mytilelist in this example) containing multiple, ">&2 
	echo "white space delimited tile IDs can be given in place of tileid ">&2
	echo "to download for multiple tiles." >&2
     	exit 1
fi

set -u
set -o pipefail

### Determine wheter a tile ID is given or a file of tile IDs is given. Validate the parameters.
case $tilePara in 
  [0-6][0-9][A-Z][A-Z][A-Z]) all_tileids=$tilePara;;
  *)   	if [ -f $tilePara ]
	then
		all_tileids=$(cat $tilePara)
	else
		echo "Parameter $tilePara is neither a tile ID nor a text file of tile IDs. exit" >&2
		exit 1
	fi
	;;
esac


### wget or curl?
which wget >/dev/null 2>&1
WGET=$? 
which curl >/dev/null 2>&1 
CURL=$? 

if [ $WGET -ne 0 ] && [ $CURL -ne 0 ]
then
	echo "Neither wget nor curl is available on your system. Exit" >&2
	exit 1
fi


datapath=https://hls.gsfc.nasa.gov/data/v1.4

################################################################################
### A function parsing the html for the items to download. 
#   Maybe some feature of the powerful wget unknown to me can make this function entirely unnecessary.
function parse_html
{
 # <tr><td valign="top"><img src="/icons/folder.gif" alt="[DIR]"></td><td><a href="2017/">2017/</a></td> ......
        local html=$1
        items=$( awk '	BEGIN {n = 1000000000 }  # an arbitrary big number
			$0 ~ /Parent Directory/ {n = NR}
			NR > n && $0 ~ /href/ { item = ""
						pos = index($0, "href")
						for (i = 0; ; i++) {
							c = substr($0, pos+6+i, 1)
							if (c == "/" || c == "\"")
								break
							else
								item = item c
						}
						print item
					      }
		     ' $html )
        echo $items
}

html=$OUTDIR/tmp.index.html
producthtml=$OUTDIR/tmp_all_hdf.html
outfile=
trap 'rm -f $html $producthtml' 0
trap ' 	case "$outfile" in
	  *hdf) base=$(basename $outfile);;
          *hdr) base=$(basename $outfile .hdr);;
	esac

      	dir=$(dirname $outfile)
      	rm -f $dir/${base}*; 

	exit 1 ' 1 2 15 	# Delete partially downloaded files on interruption. Delete both hdf and hdr

### The HLS data types as indicated by the subdirectory names in $datapath. 
#   They should be L30 and S30 for now.

if [ $WGET -eq 0 ] 
then
	wgetoption="--no-check-certificate" 
	wget  $datapath/ $wgetoption -q -O $html     
	if [ $? -ne 0 ]
	then
		echo "Is this a valid URL?  $datapath" >&2
		exit 1
	fi
	
	for hlstype in $( parse_html $html)
	do
		for tileid in $all_tileids
		do
			subdir=$(echo $tileid | awk '{print substr($1,1,2) "/" substr($1,3,1) "/" substr($1,4,1) "/" substr($1,5,1) }')
		
			if [ -n "$YEAR" ]
			then
				allyears=$YEAR
			else
				wget  $datapath/$hlstype $wgetoption -q -O $html
				allyears=$(parse_html $html)
			fi

			for year in $allyears
			do
				outdir=$OUTDIR
				mkdir -p $outdir 
				# No "-q" here; be talkative

				# Sep 20, 2018: The resursive option turned out not be as appealing because unrelated files may show up on the
				# screen before wget screens them, although they are not downloaded. Causing confusion.
				# wget -e robots=off --recursive --continue -nH --cut-dirs=9 -A hdf,hdr  $datapath/$hlstype/$year/$subdir/ $wgetoption  -P $outdir

				# New on Sep 20, 2018
				wget  $datapath/$hlstype/$year/$subdir/  $wgetoption -q -O $producthtml
				for product in $(parse_html $producthtml)
				do
					# Oct 10, 2018: The --continue option will bypass time stamp checking; do not use it.
					# wget -e robots=off --continue -N $datapath/$hlstype/$year/$subdir/$product  $wgetoption -P $outdir 
					#
					# The -O option also bypasses time stamp checking. So continue using -P option. Variable outfile is
					# created to delete partially downloaded files upon interruption and not used in wget.
					# wget -e robots=off -N $datapath/$hlstype/$year/$subdir/$product  $wgetoption -O $outfile
					#echo "====== Trying $product" >&2
					outfile=$outdir/$product
					wget -e robots=off -N $datapath/$hlstype/$year/$subdir/$product  $wgetoption -P $outdir -q
				done 
				
			done
		done
	done
else
	# Use curl
	CURLOPTION="--insecure --remote-time"
	curloption="$CURLOPTION"

	curl $datapath/ -s  $curloption -o $html
	if [ $? -ne 0 ]
	then
		echo "Is this a valid URL?  $datapath" >&2
		exit 1
	fi

	for hlstype in $( parse_html $html)
	do
		for tileid in $all_tileids
		do
			# Decided not to use deep sub-directory structure for output. Sep 25, 2018
			subdir=$(echo $tileid | awk '{print substr($1,1,2) "/" substr($1,3,1) "/" substr($1,4,1) "/" substr($1,5,1) }')
		
			# Years 
			if [ -n "$YEAR" ]
			then
				allyears=$YEAR
			else
				curl $datapath/$hlstype/ -s $curloption  -o $html
				allyears=$(parse_html $html)
			fi

			for year in $allyears
			do
				outdir=$OUTDIR
				mkdir -p $outdir    # Create directory, although curl can create dir by itself

				# hdf and headers. Curl can't do recursive downloading and we don't want it.
				curl $datapath/$hlstype/$year/$subdir/  -s $curloption -o $producthtml
				for product in $(parse_html $producthtml)
				do
					# Do not use the -C option which bypasses time stamp checking.
					echo "====== Trying $product" >&2
					outfile=$outdir/$product
					if [ -f $outfile  ]
					then
						curloption="$CURLOPTION -z $outfile"	# Download only if newer
					else
						curloption="$CURLOPTION"
					fi
					curl $datapath/$hlstype/$year/$subdir/$product $curloption -o $outfile
				done
			done
		done
	done
fi

exit 0
