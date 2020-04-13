SDL_STUFF := -D_REENTRANT -lSDL2 -I/usr/include/SDL2 
main:
	g++ clustering_game.cc -O4 $(SDL_STUFF) -lopenmc -o clustering_game
