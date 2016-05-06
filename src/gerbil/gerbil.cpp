
#include "../../include/gerbil/Application.h"

using namespace gerbil;

int main(int argc, char** argv) {

	//Application
	Application application;
	application.process(argc, argv);

	std::cout << "exit" << std::endl;

	return 0;
}
