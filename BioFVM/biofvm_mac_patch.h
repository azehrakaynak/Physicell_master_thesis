#ifndef BIOFVM_MAC_PATCH_H
#define BIOFVM_MAC_PATCH_H

// Fix for GCC not recognizing _Alignof on macOS headers
#if defined(__GNUC__) && !defined(__clang__)
  #ifndef _Alignof
    #define _Alignof(type) __alignof__(type)
  #endif
#endif

#endif // BIOFVM_MAC_PATCH_H
