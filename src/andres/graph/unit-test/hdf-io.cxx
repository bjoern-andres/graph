#include <stdexcept>
#include <cstdio>

#include <andres/graph/hdf5/graph.hxx>
#include <andres/graph/hdf5/digraph.hxx>
#include <andres/graph/hdf5/complete-graph.hxx>
#include <andres/graph/hdf5/grid-graph.hxx>


inline void test(bool pred) { 
    if(!pred)
        throw std::runtime_error("Test failed."); 
}
#define testThrowsException(PRED, EXC) { \
        {\
            bool threw = false; \
            try { \
                PRED; \
            } catch (EXC& e) {\
                threw = true; \
            } catch (...) {\
                throw std::runtime_error("Test threw different exception type.");\
            }\
            if(!threw) \
                throw std::runtime_error("Test did not throw exception");\
        }\
    }
#define testThrows(PRED) { \
        {\
            bool threw = false; \
            try { \
                PRED; \
            } catch (...) {\
                threw = true; \
            }\
            if(!threw) \
                throw std::runtime_error("Test did not throw exception");\
        }\
    }
    
#define testHDF5ThrowsException(PRED, EXC) { \
        {\
            bool threw = false; \
            H5E_auto2_t  oldPrinter; \
            void *oldClientData; \
            H5Eget_auto(H5E_DEFAULT, &oldPrinter, &oldClientData); \
            H5Eset_auto(H5E_DEFAULT, NULL, NULL); \
            try { \
                PRED; \
            } catch (EXC& e) {\
                threw = true; \
            } catch (...) {\
                H5Eset_auto(H5E_DEFAULT, oldPrinter, oldClientData); \
                throw std::runtime_error("Test threw different exception type.");\
            }\
            H5Eset_auto(H5E_DEFAULT, oldPrinter, oldClientData); \
            if(!threw) \
                throw std::runtime_error("Test did not throw exception");\
        }\
    }

using namespace andres::graph;

template<typename GRAPH_A, typename GRAPH_B>
void testGraphsEqual(const GRAPH_A& gA, const GRAPH_B& gB) {
    test(gA.numberOfEdges() == gB.numberOfEdges());
    test(gA.numberOfVertices() == gB.numberOfVertices());
    test(gA.multipleEdgesEnabled() == gB.multipleEdgesEnabled());
    for(std::size_t s=0;s<gA.numberOfVertices();++s) {
        for(std::size_t t=0;t<gA.numberOfVertices();++t) {
            std::pair<bool, std::size_t> qrA = gA.findEdge(s, t);
            std::pair<bool, std::size_t> qrB = gB.findEdge(s, t);
            test(qrA.first == qrB.first);
        }
    }
}

template<class T>
inline void testScalar(hid_t fileHandle, const std::string& datasetName, const T& data) {
    T readData;
    hdf5::load<T>(fileHandle, datasetName, readData);
    test(data == readData);
}

void testLowLevel(const std::string& fileName) {
    // Test File Creation
    {
        hdf5::HandleCheck<true> handleCheck;
        {
            hdf5::HandleCheck<true> handleCheck;
            hid_t fileHandle = hdf5::createFile(fileName);
            hdf5::closeFile(fileHandle);
        }
        {
            hdf5::HandleCheck<true> handleCheck;
            hid_t fileHandle = hdf5::openFile(fileName);
            hdf5::closeFile(fileHandle);
        }
    }
    // Test file group creation
    {
        // Group Creation
        {
            hdf5::HandleCheck<true> handleCheck;
            hid_t fileHandle = hdf5::createFile(fileName);
            hid_t groupHandle = hdf5::createGroup(fileHandle, "group");
            hdf5::closeGroup(groupHandle);
            hdf5::closeFile(fileHandle);
        }
        // Group Create/Open
        {
            hdf5::HandleCheck<true> handleCheck;
            hid_t fileHandle = hdf5::openFile(fileName, hdf5::FileAccessMode::READ_ONLY);
            hid_t groupHandle = hdf5::openGroup(fileHandle, "group");
            hdf5::closeGroup(groupHandle);
            hdf5::closeFile(fileHandle);
        }
        // File Access
        {
            hdf5::HandleCheck<true> handleCheck;
            hid_t fileHandle = hdf5::openFile(fileName, hdf5::FileAccessMode::READ_ONLY);
            testHDF5ThrowsException(hdf5::createGroup(fileHandle, "group-new"), std::runtime_error);
            hdf5::closeFile(fileHandle);
        }
        // Group OpenOrCreate
        {
            hdf5::HandleCheck<true> handleCheck;
            hid_t fileHandle = hdf5::createFile(fileName);
            {
                hdf5::HandleCheck<true> handleCheck;
                testHDF5ThrowsException(hdf5::openGroup(fileHandle, "group"), std::runtime_error);
                testHDF5ThrowsException(hdf5::openGroup(fileHandle, "group", false), std::runtime_error);
            }
            hdf5::closeFile(fileHandle);
        }
    }
    // Test scalars
    {
        const char vChar = (1<<8)-1;
        const short vShort = 2;
        const long vLong = 3;
        const long long vLongLong = 4;
        const float vFloat = 2.0/3;
        const double vDouble = 1.0/3;
        
        {
            hdf5::HandleCheck<true> handleCheck;
            hid_t fileHandle = hdf5::createFile(fileName);
            // write scalars
            hdf5::save<unsigned char>(fileHandle, "scalar-uchar", vChar);
            hdf5::save<unsigned short>(fileHandle, "scalar-ushort", vShort);
            hdf5::save<unsigned long>(fileHandle, "scalar-ulong", vLong);
            hdf5::save<unsigned long long>(fileHandle, "scalar-ulonglong", vLongLong);
            hdf5::save<signed char>(fileHandle, "scalar-char", -vChar);
            hdf5::save<signed short>(fileHandle, "scalar-short", -vShort);
            hdf5::save<signed long>(fileHandle, "scalar-long", -vLong);
            hdf5::save<signed long long>(fileHandle, "scalar-longlong", -vLongLong);
            hdf5::save<float>(fileHandle, "scalar-float", vFloat);
            hdf5::save<double>(fileHandle, "scalar-double", vDouble);
            hdf5::closeFile(fileHandle);
        }
        {
            hdf5::HandleCheck<true> handleCheck;
            hid_t fileHandle = hdf5::openFile(fileName, hdf5::FileAccessMode::READ_ONLY);
            testScalar<unsigned char>(fileHandle, "scalar-uchar", vChar);
            testScalar<unsigned short>(fileHandle, "scalar-ushort", vShort);
            testScalar<unsigned long>(fileHandle, "scalar-ulong", vLong);
            testScalar<unsigned long long>(fileHandle, "scalar-ulonglong", vLongLong);
            testScalar<signed char>(fileHandle, "scalar-char", -vChar);
            testScalar<signed short>(fileHandle, "scalar-short", -vShort);
            testScalar<signed long>(fileHandle, "scalar-long", -vLong);
            testScalar<signed long long>(fileHandle, "scalar-longlong", -vLongLong);
            testScalar<float>(fileHandle, "scalar-float", vFloat);
            testScalar<double>(fileHandle, "scalar-double", vDouble);
            // Test exceptions
            {
                signed char readChar;
                unsigned short readUShort;
                float readFloat;
                testThrowsException(hdf5::load<signed char>(fileHandle, "scalar-uchar", readChar), std::runtime_error);
                testThrowsException(hdf5::load<unsigned short>(fileHandle, "scalar-char", readUShort), std::runtime_error);
                testThrowsException(hdf5::load<float>(fileHandle, "scalar-char", readFloat), std::runtime_error);
            }
            hdf5::closeFile(fileHandle);
        }
        
    }
}

template<class GRAPH>
void saveGraph(const std::string& fileName, const std::string& datasetName, const GRAPH& g) {
    hid_t fileHandle = hdf5::createFile(fileName);
    try{
        hdf5::save(fileHandle, datasetName, g);
    } catch (std::exception& e) {
        hdf5::closeFile(fileHandle);
        throw std::runtime_error("test-hdf: Saving dataset failed: " + std::string(e.what()));
    }
    hdf5::closeFile(fileHandle);
}

template<class GRAPH>
void loadGraph(const std::string& fileName, const std::string& datasetName, GRAPH& g) {
    hid_t fileHandle = hdf5::openFile(fileName);
    try{
        hdf5::load(fileHandle, datasetName, g);
    } catch (std::exception& e) {
        hdf5::closeFile(fileHandle);
        throw std::runtime_error("test-hdf: Reading dataset failed: " + std::string(e.what()));
    }
    hdf5::closeFile(fileHandle);
}

void testSaveLoadGraph(const std::string& fileName) {
    {
        std::size_t vecI[] = {1, 1, 1, 1, 2, 2, 2, 3, 3, 4, 5, 5};
        std::size_t vecJ[] = {1, 2, 3, 4, 3, 5, 6, 4, 5, 5, 6, 7};
        
        std::size_t numberOfVertices = 10;
        
        const std::size_t numberOfEdges = sizeof(vecI)/sizeof(vecI[0]);
        
        Graph<> g(numberOfVertices);
        for(std::size_t e=0;e<numberOfEdges;++e) {
            g.insertEdge(vecI[e], vecJ[e]);
        }
        Graph<> gLoaded;
        saveGraph(fileName, "graph", g);
        loadGraph(fileName, "graph", gLoaded);
        testGraphsEqual(g, gLoaded);
    }
    // Test multiple Edges
    {
        std::size_t vecI[] = {0, 1, 1, 1, 1, 1, 1, 4, 5, 6};
        std::size_t vecJ[] = {1, 2, 3, 4, 1, 1, 1, 5, 4, 2};
        
        std::size_t numberOfVertices = 10;
        
        const std::size_t numberOfEdges = sizeof(vecI)/sizeof(vecI[0]);
        
        Graph<> g(numberOfVertices);
        g.multipleEdgesEnabled() = true;
        for(std::size_t e=0;e<numberOfEdges;++e) {
            g.insertEdge(vecI[e], vecJ[e]);
        }
        Graph<> gLoaded;
        saveGraph(fileName, "graph", g);
        loadGraph(fileName, "graph", gLoaded);
        testGraphsEqual(g, gLoaded);
    }
}

void testSaveLoadDigraph(const std::string& fileName) {
    {
        std::size_t vecI[] = {1, 1, 1, 2, 2, 2, 3, 3, 4, 5, 5, 2, 3, 5, 7};
        std::size_t vecJ[] = {2, 3, 4, 3, 5, 6, 4, 5, 5, 6, 7, 1, 1, 2, 3};
        
        std::size_t numberOfVertices = 10;
        
        const std::size_t numberOfEdges = sizeof(vecI)/sizeof(vecI[0]);
        
        Digraph<> g(numberOfVertices);
        for(std::size_t e=0;e<numberOfEdges;++e) {
            g.insertEdge(vecI[e], vecJ[e]);
        }
        Digraph<> gLoaded;
        saveGraph(fileName, "graph", g);
        loadGraph(fileName, "graph", gLoaded);
        testGraphsEqual(g, gLoaded);
    }
    // Test multiple edges
    {
        std::size_t vecI[] = {1, 1, 1, 2, 2, 2, 3, 3, 4, 5, 5, 2, 3, 5, 7, 1, 1, 3, 3, 3};
        std::size_t vecJ[] = {2, 3, 4, 3, 5, 6, 4, 5, 5, 6, 7, 1, 1, 2, 3, 2, 2, 1, 1, 1};
        
        std::size_t numberOfVertices = 10;
        
        const std::size_t numberOfEdges = sizeof(vecI)/sizeof(vecI[0]);
        
        Digraph<> g(numberOfVertices);
        g.multipleEdgesEnabled() = true;
        for(std::size_t e=0;e<numberOfEdges;++e) {
            g.insertEdge(vecI[e], vecJ[e]);
        }
        Digraph<> gLoaded;
        saveGraph(fileName, "graph", g);
        loadGraph(fileName, "graph", gLoaded);
        testGraphsEqual(g, gLoaded);
    }
}

void testSaveLoadGridGraph(const std::string& fileName) {
    // Test simple IO
    {
        GridGraph<3> g({3, 4, 2});
        GridGraph<3> gLoaded;
        saveGraph(fileName, "graph", g);
        loadGraph(fileName, "graph", gLoaded);
        testGraphsEqual(g, gLoaded);
    }
    // Test Errors
    {
        GridGraph<3> g({3, 4, 2});
        GridGraph<2> gLoaded;
        saveGraph(fileName, "graph", g);
        testThrowsException(loadGraph(fileName, "graph", gLoaded), std::runtime_error);
    }
}

void testSaveLoadCompleteGraph(const std::string& fileName) {
    typedef CompleteGraph<> CompleteGraph;
    CompleteGraph g(20);
    CompleteGraph gLoaded;
    saveGraph(fileName, "graph", g);
    loadGraph(fileName, "graph", gLoaded);
    testGraphsEqual(g, gLoaded);
}

void testTypeDetection(const std::string& fileName) {
    {
        saveGraph(fileName, "graph", Graph<>());
        { Digraph<> g; testThrowsException(loadGraph(fileName, "graph", g), std::runtime_error); }
        { CompleteGraph<> g; testThrowsException(loadGraph(fileName, "graph", g), std::runtime_error); }
        { GridGraph<> g; testThrowsException(loadGraph(fileName, "graph", g), std::runtime_error); }
    }
    {
        saveGraph(fileName, "graph", Digraph<>());
        { Graph<> g; testThrowsException(loadGraph(fileName, "graph", g), std::runtime_error); }
        { CompleteGraph<> g; testThrowsException(loadGraph(fileName, "graph", g), std::runtime_error); }
        { GridGraph<> g; testThrowsException(loadGraph(fileName, "graph", g), std::runtime_error); }
    }
    {
        saveGraph(fileName, "graph", CompleteGraph<>(10));
        { Graph<> g; testThrowsException(loadGraph(fileName, "graph", g), std::runtime_error); }
        { Digraph<> g; testThrowsException(loadGraph(fileName, "graph", g), std::runtime_error); }
        { GridGraph<> g; testThrowsException(loadGraph(fileName, "graph", g), std::runtime_error); }
    }
    {
        saveGraph(fileName, "graph", GridGraph<3>({2, 3, 4}));
        { Graph<> g; testThrowsException(loadGraph(fileName, "graph", g), std::runtime_error); }
        { Digraph<> g; testThrowsException(loadGraph(fileName, "graph", g), std::runtime_error); }
        { CompleteGraph<> g; testThrowsException(loadGraph(fileName, "graph", g), std::runtime_error); }
    }
}

int main() {
    // const std::string fileName = tmpnam(nullptr);
    const std::string fileName = "temp.h5";
    testLowLevel(fileName);
    testSaveLoadGraph(fileName);
    testSaveLoadDigraph(fileName);
    testSaveLoadGridGraph(fileName);
    testSaveLoadCompleteGraph(fileName);
    testTypeDetection(fileName);
    return 0;
}
