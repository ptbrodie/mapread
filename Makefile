mapread: mapsrc/mapread.c mapsrc/mapread.h
	gcc -g -o mapread mapsrc/mapuser.c mapsrc/mapread.c iosrc/fileio.c alignsrc/align.c sfxsrc/suffix.c

clean:
	/bin/rm -rf mapread mapread.dSYM
