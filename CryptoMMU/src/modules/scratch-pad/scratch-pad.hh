/*
  scratch-pad.h - sample code for a Simics device
*/

#ifndef SCRATCH_PAD_H
#define SCRATCH_PAD_H

#include <stdint.h>
#include <vector>
#include <iostream>
class ScratchPad
 {
public:
  std::vector<uint8_t> data;
  std::vector<bool> valid;
  
  void CheckSize(uint64_t size)
  {
    //std::cout << "scratch pad size" << size << std::endl;
    while (data.size() <= size) {
      //std::cout << "size of vector" << data.size() << std::endl;
      data.push_back(0xcd);
      valid.push_back(false);
    }
  }

  void Clear()
  {
    data = std::vector<uint8_t>();
    valid = std::vector<bool>();
  }
};

//size is 8 bytes because only one pointer
typedef struct ScratchPadHandle_t {
  ScratchPad* sp;
} ScratchPadHandle;

//size is 24 bytes because only three pointer
typedef struct scratch_pad_interface {
  void (*read)(ScratchPadHandle* obj, uint64_t address, void* dataRd,
               unsigned int size);
  void (*write)(ScratchPadHandle* obj, uint64_t address, const void* dataWr,
                unsigned int size);
  void (*clear)(ScratchPadHandle* obj);
} scratch_pad_interface_t;

ScratchPadHandle* CreateNewScratchPad();
int DeleteScratchPad(ScratchPadHandle* obj);
scratch_pad_interface_t* CreateScratchPadInterface();

#endif /* SCRATCH_PAD_H */
