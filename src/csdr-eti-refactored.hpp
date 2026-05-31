// Minimal buffer interface (no csdr dependency)
class BufferReader {
public:
    virtual ~BufferReader() = default;
    virtual size_t available() const = 0;
    virtual const void* getReadPointer() const = 0;
    virtual void advance(size_t bytes) = 0;
};

class BufferWriter {
public:
    virtual ~BufferWriter() = default;
    virtual size_t writeable() const = 0;
    virtual void* getWritePointer() = 0;
    virtual void advance(size_t bytes) = 0;
};

namespace Csdr::Eti {
    class EtiDecoder {
    public:
        void setBufferReader(BufferReader* reader) { this->reader = reader; }
        void setBufferWriter(BufferWriter* writer) { this->writer = writer; }

    private:
        BufferReader* reader;
        BufferWriter* writer;
    };
}
