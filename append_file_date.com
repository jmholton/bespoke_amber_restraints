#! /bin/csh -f
#
#
#

foreach file ( $* )
    if(! -e "$file") then
	echo "${file} does not exist"
	continue
    endif
    set date = `date +%m-%d-%y --reference=$file`
gotdate:
    if(-e "${file}.$date") then
        echo "${file}.$date already exists"
	diff -q $file ${file}.$date
	if($status) then
	    echo "and it is different"
	    if(! $?DATE_TIME) then
		set DATE_TIME
	        set date = `date +%m-%d-%y_%H:%M:%S --reference=$file`
	        goto gotdate
	    endif
	else
	    echo "but identical to $file"
	endif
    else
        echo "cp -p $file ${file}.$date"
        cp -p $file ${file}.$date
    endif
end

