//=============================================================================
#include <Inventor/Win/SoWin.h>
#include <Inventor/Win/viewers/SoWinExaminerViewer.h>
#include "MeshViewer.h"
#include "SkelatonExtractor.h"
#include "Windows.h"
#include <string>
#include <codecvt>
#include <locale>
#include "Mesh.h"
#include "Painter.h"

void loadviewer(string fname);

int main(int argc, char **argv)
{
    std::string meshName = "man0.off";
    const char *fileName = meshName.c_str();
    SkelatonExtractor mv;
    mv.open_mesh(fileName);
    //mv.draw();

    HWND window = SoWin::init(argv[0]);

	SoWinExaminerViewer * viewer = new SoWinExaminerViewer(window);

	//make a dead simple scene graph by using the Coin library, only containing a single cone under the scenegraph root
	SoSeparator * root = new SoSeparator;
	root->ref();

	Mesh* mesh = new Mesh();
	Painter* painter = new Painter();
	char *newName = new char[mv.getNewFileName().size()+1];
	strcpy(newName, mv.getNewFileName().c_str());
	mesh->loadOff("new_man0.off");
	mesh->do_edge_collapse_process();
	//mesh->writeMesh("collapsed_man0.off");
	//loadviewer("collapsed_man0.off");

	root->addChild( painter->getShapeSep(mesh) );


	viewer->setSize(SbVec2s(640, 480));
	viewer->setSceneGraph(root);
	viewer->show();

	SoWin::show(window);
	SoWin::mainLoop();
	delete viewer;
	root->unref();
    return 1;
}

void loadviewer(string fname)
{
    STARTUPINFO si;
    PROCESS_INFORMATION pi; 
    ZeroMemory(&si, sizeof(si));
    si.cb = sizeof(si);
    ZeroMemory(&pi, sizeof(pi));

    std::wstring arguments = std::wstring_convert<std::codecvt_utf8<wchar_t>>().from_bytes("0 " + fname);
    CONST wchar_t *commandLine = arguments.c_str();

    SECURITY_ATTRIBUTES secattr;
    ZeroMemory(&secattr, sizeof(secattr));
    secattr.nLength = sizeof(secattr);
    secattr.bInheritHandle = TRUE;

    HANDLE rPipe, wPipe;

    //Create pipes to write and read data
    CreatePipe(&rPipe, &wPipe, &secattr, 0);
    si.cb = sizeof(si);
    si.dwFlags = STARTF_USESTDHANDLES;
    si.hStdInput = NULL;
    si.hStdOutput = wPipe;
    si.hStdError = wPipe;

    // Start libigl viewer.
    if (!CreateProcessW(
            L"D:\\belgeler\\merve\\s02\\ceng789\\term project\\Project\\skeletonextraction\\viewer.exe", // app path
            (LPWSTR)commandLine,                                                                         // Command line
            NULL,                                                                                        // Process handle not inheritable
            NULL,                                                                                        // Thread handle not inheritable
            FALSE,                                                                                       // Set handle inheritance to FALSE
            0,                                                                                           // No creation flags
            NULL,                                                                                        // Use parent's environment block
            NULL,                                                                                        // Use parent's starting directory
            &si,                                                                                         // Pointer to STARTUPINFO structure
            &pi)                                                                                         // Pointer to PROCESS_INFORMATION structure
    )
    {
        printf("CreateProcess failed (%d).\n", GetLastError());
        throw std::exception("Could not create viwer process");
    }
    else
    {
        std::cout << "[          ] Successfully launched child process" << std::endl;
    }
}