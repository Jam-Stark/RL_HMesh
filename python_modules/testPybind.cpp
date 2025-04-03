//test if pybind11 lauched
#include <pybind11/embed.h>
#include <iostream>

namespace py = pybind11;

using namespace std;

int main(int argc, char *argv[]) {
    py::scoped_interpreter guard{};
    
    try {
        auto sys = py::module_::import("sys");
        sys.attr("path").attr("append")(".");
        sys.attr("path").attr("append")(R"(F:\RL_HMesh\.venv\Lib\site-packages)");  // 关键路径添加
        
        // 调试输出
        std::cout << "Using Python: " << sys.attr("executable").cast<std::string>() << "\n";
        std::cout << "Python path: \n";
        for (auto p : sys.attr("path")) {
            std::cout << "  " << p.cast<std::string>() << "\n";
        }

        auto test_module = py::module_::import("test");
        test_module.attr("train")();
        
    } catch (const py::error_already_set &e) {
        std::cerr << "Python Error:\n" << e.what() << std::endl;
    } catch (const std::exception &e) {
        std::cerr << "C++ Error: " << e.what() << std::endl;
    }

    return 0;
}


