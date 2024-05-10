#Write inputRun.xml file for a cubic domain of NxNxN blocks (each block is 8x8x8 sites).
# Run as bash writeInput.sh N

cen=$(( (8*$1 +1) / 2))
end=$(( 8*$1 ))

awk -v A=$cen -v B=$end '{gsub(/CEN/,A);sub(/END/,B)}1' inputBase.xml > inputRun.xml

