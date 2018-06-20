# tcsh host and Makefile autocompletion (tab expansions)

# shamelessly stolen from 
# Debian GNU/Linux
# /usr/share/doc/tcsh/examples/complete.gz
#
# See that file for copyright notices (and lots more completions).


# intro sanity check
onintr -
if (! $?prompt) goto end

# get some global variables and version info
if ($?tcsh) then
    if ($tcsh != 1) then
   	set rev=$tcsh:r
	set rel=$rev:e
	set pat=$tcsh:e
	set rev=$rev:r
    endif
    if ($rev > 5 && $rel > 1) then
	set _complete=1
    endif
    unset rev rel pat
endif

# now for the actual autocompletion
if ($?_complete) then
    set noglob
    #
    # tab expansion for hosts (used for ssh etc.)
    #
    #
    # step 1: scan a couple of files where host information is usually stored
    if ( ! $?hosts ) set hosts
    foreach f ("$HOME/.hosts" /usr/local/etc/csh.hosts "$HOME/.rhosts" /etc/hosts.equiv)
        if ( -r "$f" ) then
	    set hosts = ($hosts `grep -v "+" "$f" | grep -E -v "^#" | tr -s " " "	" | cut -f 1`)
	endif
    end
    if ( -r "$HOME/.netrc" ) then
	set f=`awk '/machine/ { print $2 }' < "$HOME/.netrc"` >& /dev/null
	set hosts=($hosts $f)
    endif
    if ( -r "$HOME/.ssh/known_hosts" ) then
	set f=`cat "$HOME/.ssh/known_hosts" | cut -f 1 -d \ ` >& /dev/null
	set f=`cat "$HOME/.ssh/known_hosts" | cut -f 1 -d \ | sed -e 's/,/ /g'` >& /dev/null
	set hosts=($hosts $f)
    endif
    unset f
    # step 2: add a bunch of defaults (TODO: Add more)
    if ( ! $?hosts ) then
	set hosts=(teslaspule.mathematik.tu-dortmund.de \
		   kittyhawk.mathematik.tu-dortmund.de)
    endif
    # step 3: enable tab expansions
    complete ywho  	n/*/\$hosts/	# argument from list in $hosts
    complete rsh	p/1/\$hosts/ c/-/"(l n)"/   n/-l/u/ N/-l/c/ n/-/c/ p/2/c/ p/*/f/
    complete ssh	p/1/\$hosts/ c/-/"(l n)"/   n/-l/u/ N/-l/c/ n/-/c/ p/2/c/ p/*/f/
    complete xrsh	p/1/\$hosts/ c/-/"(l 8 e)"/ n/-l/u/ N/-l/c/ n/-/c/ p/2/c/ p/*/f/
    complete rlogin 	p/1/\$hosts/ c/-/"(l 8 e)"/ n/-l/u/
    complete telnet 	p/1/\$hosts/ p/2/x:'<port>'/ n/*/n/
    complete finger	c/*@/\$hosts/ n/*/u/@ 
    complete ping	p/1/\$hosts/
    complete traceroute	p/1/\$hosts/
    complete nslookup   p/1/x:'<host>'/ p/2/\$hosts/
    #
    #
    # autocompletion for make and gmake
    #
    #
    complete gmake	'c/{--directory=,--include-dir=}/d/' \
			'c/{--assume-new,--assume-old,--makefile,--new-file,--what-if,--file}/f/' \
			'c/--/(assume-new= assume-old= debug directory= \
			dry-run environment-overrides file= help \
			ignore-errors include-dir= jobs[=N] just-print \
			keep-going load-average[=N] makefile= max-load[=N] \
			new-file= no-builtin-rules no-keep-going \
			no-print-directory old-file= print-data-base \
			print-directory question quiet recon silent stop \
			touch version warn-undefined-variables what-if=)/' \
			'n@*@`cat -s GNUMakefile Makefile makefile |& sed -n -e "/No such file/d" -e "s/^\([A-Za-z0-9-]*\):.*/\1/p"`@' \
			'n/=/f/' 'n/-f/f/'
    complete make	'c/{--directory=,--include-dir=}/d/' \
			'c/{--assume-new,--assume-old,--makefile,--new-file,--what-if,--file}/f/' \
			'c/--/(assume-new= assume-old= debug directory= \
			dry-run environment-overrides file= help \
			ignore-errors include-dir= jobs[=N] just-print \
			keep-going load-average[=N] makefile= max-load[=N] \
			new-file= no-builtin-rules no-keep-going \
			no-print-directory old-file= print-data-base \
			print-directory question quiet recon silent stop \
			touch version warn-undefined-variables what-if=)/' \
			'n@*@`cat -s GNUMakefile Makefile makefile |& sed -n -e "/No such file/d" -e "s/^\([A-Za-z0-9-]*\):.*/\1/p"`@' \
			'n/=/f/' 'n/-f/f/'


    # done, clean up
    unset noglob
    unset _complete
    unset traditional_complete
endif

end:
	onintr
