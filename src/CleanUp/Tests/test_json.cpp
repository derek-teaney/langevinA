#include "json/json.h"
#include <cstdio>
#include <fstream>
#include <iostream>

int main(int argc, char **argv) {
  printf("Hello World\n");
  Json::Value parser;
  std::ifstream ifs("test_json.json");
  ifs >> parser;

  int NX = parser["NX"].asInt();
  int NY = parser.get("NY", NX).asInt();
  int NZ = parser["NZ"].asInt();

  std::cout << "NX=" << NX << std::endl;
  std::cout << "NY=" << NY << std::endl;
  std::cout << "NZ=" << NZ << std::endl;

  std::cout << parser << std::endl;
}
