Windows
==================================

2.Configurate CMake on Windows
	We could build the C++ compiler with MinGW.Here is the link
	
	https://sourceforge.net/projects/mingw-w64/
	
	Running the setup file and set the architecture to "x86_64", also remember the installation directory.	

	Then go to the installation folder, add the ``<mingw installation path>\i686-8.1.0-posix-dwarf-rt_v6-rev0/mingw32/include`` to the windows environment path. By now, MinGW should have been configurated successfully. 
	
	MinGW has a package manager which would be convenient to install the libraries and packages. We could setup OpenMP with it.
	
	https://sourceforge.net/projects/mingw/
	
	The installation process is quite straight forward. After installation, you could be able to use "MinGW Installation Manager". Some packages and libraries could be installed directly with this software
	
.. image:: windows3.PNG


3. cmake and compile
	Open the cmd and go to the directroy of the "CmakeLists.txt". Then
	
	``cmake -G "Unix Makefiles"``
	
	``make``
	
	The executable file will be generated to the "dist/Debug/GNU-Linux/nutrient_app" folder.
	
