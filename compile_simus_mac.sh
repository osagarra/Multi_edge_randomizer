clear
clang -c ./includes/stat_funcs.c -O2 -lm -lgsl -lgslcblas -Wall -g
clang -c ./includes/ula_funcs.c -O2 -lm -lgsl -lgslcblas -Wall -g
clang -c ./includes/graph_funcs.c -O2 -lm -lgsl -lgslcblas -Wall -g
clang -c ./includes/w_graph_funcs.c -O2 -lm -lgsl -lgslcblas -Wall -g
clang -c ./includes/transport_models/null_models.c -O2 -lm -lgsl -lgslcblas -Wall -g

clang ./main_simus_cust.c ula_funcs.o stat_funcs.o graph_funcs.o w_graph_funcs.o null_models.o -O2 -lm -lgsl -lgslcblas -Wall -o simus_cust -g
