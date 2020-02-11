## sort through just the fourth line of each read (fourth line is quality line)
awk '{if(NR%4==0){print $0}}' SRR2859223_1.fastq > /scratch/cat7ep/simCline/biosampleresults/SRR2859223_1out.txt

## Write an array of the ASCII scale of the quality markers
## !"#$%&'()*+,-./0123456789:;<=>?@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~
declare -a quality
quality=(\! \" \# \$ \% \& \' \( \) \* \+ \, \- \. \/ 0 1 2 3 4 5 6 7 8 9 \: \; \< \= \> \? \@ A B C D E F G H I J K L M N O P Q R S T U V W X Y Z \[ \\ \] \^ \_ \` a b c d e f g h i j k l m n o p q r s t u v w x y z \{ \| \} \~ )

## Search through each row of quality, get the range
if((grep -c ${quality[0-25]} SRR2859223_1out.txt)>0){
  echo "Sanger or Illumina 1.8+"
}
## match the range to certain read (with confidence interval?)
