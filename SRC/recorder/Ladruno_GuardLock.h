// LADRUNO-HEADER-START
// ==========================================================================
//
//   ▄█          ▄████████ ████████▄     ▄████████ ███    █▄  ███▄▄▄▄    ▄██████▄
//  ███         ███    ███ ███   ▀███   ███    ███ ███    ███ ███▀▀▀██▄ ███    ███
//  ███         ███    ███ ███    ███   ███    ███ ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███  ▄███▄▄▄▄██▀ ███    ███ ███   ███ ███    ███
//  ███       ▀███████████ ███    ███ ▀▀███▀▀▀▀▀   ███    ███ ███   ███ ███    ███
//  ███         ███    ███ ███    ███ ▀███████████ ███    ███ ███   ███ ███    ███
//  ███▌    ▄   ███    ███ ███   ▄███   ███    ███ ███    ███ ███   ███ ███    ███
//  █████▄▄██   ███    █▀  ████████▀    ███    ███ ████████▀   ▀█   █▀   ▀██████▀
//  ▀                                   ███    ███
//
//  Ladruno — a research fork of OpenSees
//  Created by:  Nicolas Mora Bowen  ·  Patricio Palacios  ·  José Abell  ·  Guppi
//
// Header auto-stamped by Ladruno_scripts/stamp_headers.py (art: banner_ASCII.txt).
// Do not hand-edit between the markers; edit the script/art and re-run instead.
// ==========================================================================
// LADRUNO-HEADER-END

#ifndef Ladruno_GuardLock_h
#define Ladruno_GuardLock_h

/*************************************************************************************

Ladruno_GuardLock.h

The guard lock of every SWMR writer of a `.ladruno` file (WP-233, #1082; shared with the
MPM particle writer by WP-232h, #1085). Moved here VERBATIM from LadrunoRecorder.cpp.

A reader with HDF5 file locking on (h5py's default, HDFView, ...) holds a shared lock
while it opens the file, and with HDF5's own locking that breaks the SWMR writer
(measured): H5Fstart_swmr_write and a layout step's reopen must lock and are refused
(Linux: EAGAIN), and on Windows, where byte-range locks are mandatory, the writer's very
writes fail ("Permission denied"). So HDF5's locking is off for the file
(Ladruno_Hdf5.h, h5::file::accessList) and the writer holds its own EXCLUSIVE lock, from
right after H5Fcreate until after the final close, through a separate descriptor:
  * Windows: LockFileEx on one byte at offset 2^62, which every reader's whole-file
    lock overlaps; the writer's own writes never reach it;
  * POSIX: flock(LOCK_EX), which every reader's flock(LOCK_SH) conflicts with; flock is
    advisory, so it never blocks the writer's writes.
A locking reader is then refused ("unable to lock file") for as long as the run is live,
as before WP-233, when the writer held HDF5's exclusive lock all run. A reader with
locking OFF (as the live monitor requires, 08_analysis_monitor.md) takes no lock and
reads the live file through SWMR. The OS releases the guard if the process is killed.

acquire() returns 0 when the file cannot be opened or a reader holds a shared lock for
the whole ~60 s retry; the caller warns once and carries on unguarded.

**************************************************************************************/

#include <algorithm>
#include <chrono>
#include <cstdint>
#include <cstring>
#include <string>
#include <thread>

#if defined(_WIN32)
#ifndef NOMINMAX
#define NOMINMAX
#endif
#ifndef WIN32_LEAN_AND_MEAN
#define WIN32_LEAN_AND_MEAN
#endif
#ifndef NOGDI
#define NOGDI
#endif
#include <windows.h>
#else
#include <fcntl.h>      // flock on a separate descriptor
#include <sys/file.h>
#include <unistd.h>
#endif

namespace ladruno {
namespace guardlock {

inline void* acquire(const std::string& path)
{
#if defined(_WIN32)
	// A locking reader may hold a shared lock this instant: retry for up to ~60 s.
	for (int attempt = 0, pause_ms = 2; attempt < 130; ++attempt,
	     pause_ms = std::min(pause_ms * 2, 500)) {
		HANDLE h = CreateFileA(path.c_str(), GENERIC_READ,
		                       FILE_SHARE_READ | FILE_SHARE_WRITE | FILE_SHARE_DELETE,
		                       NULL, OPEN_EXISTING, FILE_ATTRIBUTE_NORMAL, NULL);
		if (h == INVALID_HANDLE_VALUE)
			return 0;
		OVERLAPPED ov;
		std::memset(&ov, 0, sizeof(ov));
		ov.Offset = 0;
		ov.OffsetHigh = 0x40000000u;   // byte 2^62: far past any data
		if (LockFileEx(h, LOCKFILE_EXCLUSIVE_LOCK | LOCKFILE_FAIL_IMMEDIATELY, 0, 1, 0, &ov))
			return (void*)h;
		CloseHandle(h);   // a reader holds a shared lock this instant: retry
		std::this_thread::sleep_for(std::chrono::milliseconds(pause_ms));
	}
	return 0;
#else
	for (int attempt = 0, pause_ms = 2; attempt < 130; ++attempt,
	     pause_ms = std::min(pause_ms * 2, 500)) {
		const int fd = ::open(path.c_str(), O_RDONLY);
		if (fd < 0)
			return 0;
		if (::flock(fd, LOCK_EX | LOCK_NB) == 0)
			return (void*)(intptr_t)(fd + 1);   // + 1: fd 0 is valid, 0 means none
		::close(fd);      // a reader holds a shared lock this instant: retry
		std::this_thread::sleep_for(std::chrono::milliseconds(pause_ms));
	}
	return 0;
#endif
}

inline void release(void*& guard)
{
#if defined(_WIN32)
	if (guard) {
		OVERLAPPED ov;
		std::memset(&ov, 0, sizeof(ov));
		ov.OffsetHigh = 0x40000000u;
		UnlockFileEx((HANDLE)guard, 0, 1, 0, &ov);
		CloseHandle((HANDLE)guard);
	}
#else
	if (guard) {
		const int fd = (int)(intptr_t)guard - 1;
		::flock(fd, LOCK_UN);
		::close(fd);
	}
#endif
	guard = 0;
}

}  // namespace guardlock
}  // namespace ladruno

#endif // Ladruno_GuardLock_h
