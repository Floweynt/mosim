#include "archive.h"
#include "assets.h"
#include "resources.h"
#include <archive_entry.h>
#include <sys/stat.h>

resource_manager::resource_manager()
{
    archive* archive = archive_read_new();
    archive_read_support_format_tar(archive);
    auto r = archive_read_open_memory(archive, assets_tar, assets_tar_len);
    archive_entry* entry = archive_entry_new();

    while (archive_read_next_header(archive, &entry) == ARCHIVE_OK)
    {
        if(S_ISREG(archive_entry_filetype(entry)))
        {
            auto size = archive_entry_size(entry);
            std::vector<uint8_t> buf;
            buf.resize(size);
        
            archive_read_data(archive, (void*) buf.data(), size);
            resources[archive_entry_pathname(entry)] = std::move(buf); 
        }
    }

    archive_read_free(archive);
}

