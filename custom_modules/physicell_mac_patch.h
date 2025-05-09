#ifndef PHYSICELL_MAC_PATCH_H
#define PHYSICELL_MAC_PATCH_H

// Fix for GCC not recognizing _Alignof on macOS
#if defined(__GNUC__) && !defined(__clang__)
#ifndef _Alignof
#define _Alignof(type) __alignof__(type)
#endif
#endif

#endif // PHYSICELL_MAC_PATCH_H
