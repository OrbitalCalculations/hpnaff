//
//  hpnaff-BridgingHeader.h
//  hpnaff
//
//  Created by Heiko Pälike on 24/05/2018.
//  Copyright © 2018 Heiko Pälike. All rights reserved.
//

#ifndef hpnaff_BridgingHeader_h
#define hpnaff_BridgingHeader_h
//#import "Vendor/Surge/Sources/Surge/HP/MTRandom.h"
//#import "Vendor/Rainbow/Sources/Rainbow.h"
#endif /* hpnaff_BridgingHeader_h */

/*
#import <mach-o/dyld.h>
#import <mach-o/loader.h>

static NSUUID *ExecutableUUID(void)
{
    const struct mach_header *executableHeader = NULL;
    for (uint32_t i = 0; i < _dyld_image_count(); i++)
    {
        const struct mach_header *header = _dyld_get_image_header(i);
        if (header->filetype == MH_EXECUTE)
        {
            executableHeader = header;
            break;
        }
    }

    if (!executableHeader)
        return nil;

    BOOL is64bit = executableHeader->magic == MH_MAGIC_64 || executableHeader->magic == MH_CIGAM_64;
    uintptr_t cursor = (uintptr_t)executableHeader + (is64bit ? sizeof(struct mach_header_64) : sizeof(struct mach_header));
    const struct segment_command *segmentCommand = NULL;
    for (uint32_t i = 0; i < executableHeader->ncmds; i++, cursor += segmentCommand->cmdsize)
    {
        segmentCommand = (struct segment_command *)cursor;
        if (segmentCommand->cmd == LC_UUID)
        {
            const struct uuid_command *uuidCommand = (const struct uuid_command *)segmentCommand;
            return [[NSUUID alloc] initWithUUIDBytes:uuidCommand->uuid];
        }
    }

    return nil;
}

*/
