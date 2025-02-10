#include <fstream>
#include <iostream>
#include "ipelibrary/ipelet.h"

class MedialAxisIpelet : public ipe::Ipelet
{
public:
    virtual int ipelibVersion() const { return IPELIB_VERSION; }
    virtual bool run(int function, ipe::IpeletData *data, ipe::IpeletHelper *helper);
};

bool MedialAxisIpelet::run(int function, ipe::IpeletData *data, ipe::IpeletHelper *helper)
{
    std::ifstream infile("coordinates.txt");
    std::string polygon;

    while (getline(infile, polygon))
    {
        // Output the text from the file
        std::cout << polygon;
    }

    infile.close();
    return 0;
}

IPELET_DECLARE ipe::Ipelet *newIpelet()
{
    return new MedialAxisIpelet;
}
