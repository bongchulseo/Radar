gcc rainrate.c -m64 -Ofast -flto -march=native -funroll-loops -o rainrate -L/usr/local/lib -L../libphydro_v2 -I/usr/local/lib -I../libphydro_v2 -lphydro -lrsl -lm -lfl
