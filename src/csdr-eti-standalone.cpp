#include <iostream>
#include <fstream>
#include <cstring>
#include <cstdint>
#include <vector>
#include <getopt.h>
#include <csignal>
// #include <atomic>

extern "C" {
    #include <fftw3.h>
    #include "sdr_prstab.h"
    #include "dab_tables.h"
}

// static std::atomic<bool> g_stop_requested{false};
static volatile std::sig_atomic_t g_stop_requested = 0;
extern "C" void handle_signal(int) {
    g_stop_requested = 1;
}


// Include csdr headers
#include <csdr/complex.hpp>
#include "dab.hpp"

// Standalone decoder without csdr module system
class StandaloneEtiDecoder {
private:
    struct dab_state_t* dab = nullptr;
    uint32_t coarse_timeshift = 0;
    int32_t fine_timeshift = 0;
    int32_t coarse_freq_shift = 0;
    double fine_freq_shift = 0;
    bool force_timesync = false;

    fftwf_plan forward_plan;
    fftwf_plan backward_plan;
    fftwf_plan coarse_plan;

    std::ostream& output;
    std::vector<uint8_t> input_buffer;
    size_t read_pos = 0;

    bool sdr_demod(Csdr::complex<float>* input, struct demapped_transmission_frame_t* tf);
    uint32_t get_coarse_time_sync(Csdr::complex<float>* input);
    int32_t get_fine_time_sync(Csdr::complex<float>* input);
    int32_t get_coarse_freq_shift(Csdr::complex<float>* input);
    double get_fine_freq_corr(Csdr::complex<float>* input);

public:
    StandaloneEtiDecoder(std::ostream& out) : output(out) {
        dab = init_dab_state();

        // Set ETI callback to write directly to output
        dab->eti_callback = [this](uint8_t* eti) {
            this->output.write(reinterpret_cast<const char*>(eti), 6144);
            this->output.flush();
        };

        forward_plan = fftwf_plan_dft_1d(2048, nullptr, nullptr, FFTW_FORWARD, FFTW_ESTIMATE);
        backward_plan = fftwf_plan_dft_1d(1536, nullptr, nullptr, FFTW_BACKWARD, FFTW_ESTIMATE);
        coarse_plan = fftwf_plan_dft_1d(128, nullptr, nullptr, FFTW_BACKWARD, FFTW_ESTIMATE);
    }

    ~StandaloneEtiDecoder() {
        if (dab) delete dab;
        fftwf_destroy_plan(forward_plan);
        fftwf_destroy_plan(backward_plan);
        fftwf_destroy_plan(coarse_plan);
    }

    void addInputData(const uint8_t* data, size_t size) {
        input_buffer.insert(input_buffer.end(), data, data + size);
    }

    bool hasEnoughData() const {
        // Need 196608 complex samples = 196608 * 8 bytes
        return (input_buffer.size() - read_pos) >= (196608 * sizeof(Csdr::complex<float>));
    }

    void processFrame() {
        if (!hasEnoughData()) return;

        // Get pointer to input data
        Csdr::complex<float>* input = reinterpret_cast<Csdr::complex<float>*>(input_buffer.data() + read_pos);

        if (sdr_demod(input, &dab->tfs[dab->tfidx])) {
            dab_process_frame(dab);
        }

        // Advance input buffer
        size_t advance_samples = 196608 + coarse_timeshift + fine_timeshift;
        size_t advance_bytes = advance_samples * sizeof(Csdr::complex<float>);

        read_pos += advance_bytes;

        // Compact buffer if we've consumed more than half
        if (read_pos > input_buffer.size() / 2) {
            input_buffer.erase(input_buffer.begin(), input_buffer.begin() + read_pos);
            read_pos = 0;
        }
    }
};

// Implementation of decoder methods
bool StandaloneEtiDecoder::sdr_demod(Csdr::complex<float>* input, struct demapped_transmission_frame_t* tf) {
    coarse_timeshift = get_coarse_time_sync(input);
    force_timesync = false;
    if (coarse_timeshift) {
        std::cerr << "coarse time shift: " << coarse_timeshift << std::endl;
        return false;
    }

    if (coarse_freq_shift) {
        fine_timeshift = 0;
    } else {
        fine_timeshift = get_fine_time_sync(input);
    }

    coarse_freq_shift = get_coarse_freq_shift(input);
    if (abs(coarse_freq_shift) > 1) {
        force_timesync = true;
        return false;
    }

    fine_freq_shift = get_fine_freq_corr(input);

    /* raw symbols */
    fftwf_complex symbols[76][2048] = {0, 0};

    /* d-qpsk */
    for (int i = 0; i < 76; i++) {
        fftwf_execute_dft(forward_plan, (fftwf_complex*) &input[2656 + (2552 * i) + 504], symbols[i]);
        fftwf_complex tmp;
        for (int j = 0; j < 2048/2; j++) {
            tmp[0]     = symbols[i][j][0];
            tmp[1]     = symbols[i][j][1];
            symbols[i][j][0]    = symbols[i][j+2048/2][0];
            symbols[i][j][1]    = symbols[i][j+2048/2][1];
            symbols[i][j+2048/2][0] = tmp[0];
            symbols[i][j+2048/2][1] = tmp[1];
        }
    }

    /* symbols d-qpsk-ed */
    fftwf_complex* symbols_d = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * 2048 * 76);

    for (int j = 1; j < 76; j++) {
        for (int i = 256; i < 1793; i++) {
            symbols_d[j * 2048 + i][0] =
            ((symbols[j][i][0] * symbols[j-1][i][0])
            +(symbols[j][i][1] * symbols[j-1][i][1]))
            /(symbols[j - 1][i][0] * symbols[j - 1][i][0] + symbols[j - 1][i][1] * symbols[j - 1][i][1]);
            symbols_d[j * 2048 + i][1] =
            ((symbols[j][i][0] * symbols[j-1][i][1])
            -(symbols[j][i][1] * symbols[j-1][i][0]))
            /(symbols[j - 1][i][0] * symbols[j - 1][i][0] + symbols[j - 1][i][1] * symbols[j-1][i][1]);
        }
    }

    uint8_t* dst = tf->fic_symbols_demapped[0];

    int k, kk;
    for (int j=1; j<76; j++) {
        if (j == 4) { dst = tf->msc_symbols_demapped[0]; }
        k = 0;
        for (int i = 256; i < 1793; i++){
            if (i != 1024) {
                /* Frequency deinterleaving and QPSK demapping combined */
                kk = rev_freq_deint_tab[k++];
                dst[kk] = (symbols_d[j * 2048 + i][0] > 0) ? 0 : 1;
                dst[1536 + kk] = (symbols_d[j * 2048 + i][1] > 0) ? 1 : 0;
            }
        }
        dst += 3072;
    }

    fftwf_free(symbols_d);

    return true;
}

uint32_t StandaloneEtiDecoder::get_coarse_time_sync(Csdr::complex<float>* input) {
    int32_t tnull = 2656;
    int32_t j, k;
    std::vector<float> filt(196608 - tnull);

    float e = 0;
    float threshold = 40;
    for (k = 0; k < tnull; k += 10)
        e += fabsf(input[k].i());

    if (e < threshold && !force_timesync)
        return 0;

    memset(filt.data(), 0, sizeof(float) * (196608 - tnull) / 10);
    for (j = 0; j < (196608 - tnull) / 10; j++)
        for (k = 0; k < tnull; k += 10)
            filt[j] += fabsf(input[j * 10 + k].i());

    float minVal = 9999999;
    uint32_t minPos = 0;
    for (j = 0; j < (196608 - tnull) / 10;j++){
        if (filt[j] < minVal) {
            minVal = filt[j];
            minPos = j * 10;
        }
    }
    return minPos;
}

int32_t StandaloneEtiDecoder::get_fine_time_sync(Csdr::complex<float> *input) {
    fftwf_complex prs_received_fft[2048];
    fftwf_execute_dft(forward_plan, (fftwf_complex*) &input[2656 + 504], &prs_received_fft[0]);

    fftwf_complex prs_star[1536];
    int i;
    for (i = 0; i < 1536; i++) {
        prs_star[i][0] = Csdr::Eti::prs_static[i][0];
        prs_star[i][1] = -1 * Csdr::Eti::prs_static[i][1];
    }

    fftwf_complex prs_rec_shift[1536];
    for (i = 0; i < 1536; i++) {
        if (i < 768) {
            prs_rec_shift[i][0] = prs_received_fft[i + 1280][0];
            prs_rec_shift[i][1] = prs_received_fft[i + 1280][1];
        }
        if (i >= 768) {
            prs_rec_shift[i][0] = prs_received_fft[i - 765][0];
            prs_rec_shift[i][1] = prs_received_fft[i - 765][1];
        }
    }

    fftwf_complex convoluted_prs[1536];
    int s;
    for (s=0;s<1536;s++) {
        convoluted_prs[s][0] = prs_rec_shift[s][0] * prs_star[s][0] - prs_rec_shift[s][1] * prs_star[s][1];
        convoluted_prs[s][1] = prs_rec_shift[s][0] * prs_star[s][1] + prs_rec_shift[s][1] * prs_star[s][0];
    }

    fftwf_complex convoluted_prs_time[1536];
    fftwf_execute_dft(backward_plan, &convoluted_prs[0], &convoluted_prs_time[0]);

    int32_t maxPos=0;
    float tempVal;
    float maxVal =- 99999;
    for (i=0;i<1536;i++) {
        tempVal = sqrtf((convoluted_prs_time[i][0] * convoluted_prs_time[i][0]) + (convoluted_prs_time[i][1] * convoluted_prs_time[i][1]));
        if (tempVal > maxVal) {
            maxPos = i;
            maxVal = tempVal;
        }
    }

    if (maxPos < 1536 / 2) {
        return maxPos + 8;
    } else {
        return maxPos - 1536;
    }
}

int32_t StandaloneEtiDecoder::get_coarse_freq_shift(Csdr::complex<float> *input) {
    fftwf_complex symbols[2048] = {0};
    fftwf_execute_dft(forward_plan, (fftwf_complex*) &input[2656 + 505 + fine_timeshift], symbols);

    fftwf_complex tmp;
    for (int i = 0; i < 2048/2; i++) {
        tmp[0]     = symbols[i][0];
        tmp[1]     = symbols[i][1];
        symbols[i][0]    = symbols[i + 2048 / 2][0];
        symbols[i][1]    = symbols[i + 2048 / 2][1];
        symbols[i + 2048 / 2][0] = tmp[0];
        symbols[i + 2048 / 2][1] = tmp[1];
    }

    int len = 128;
    fftwf_complex convoluted_prs[len];
    int s;
    int freq_hub = 14;
    int k;
    float global_max = -99999;
    int global_max_pos = 0;
    for (k = -freq_hub; k <= freq_hub; k++) {
        for (s = 0; s < len; s++) {
            convoluted_prs[s][0] = Csdr::Eti::prs_static[freq_hub+s][0] * symbols[freq_hub+k+256+s][0]-
            (-1)*Csdr::Eti::prs_static[freq_hub+s][1] * symbols[freq_hub+k+256+s][1];
            convoluted_prs[s][1] = Csdr::Eti::prs_static[freq_hub+s][0] * symbols[freq_hub+k+256+s][1]+
            (-1)*Csdr::Eti::prs_static[freq_hub+s][1] * symbols[freq_hub+k+256+s][0];
        }
        fftwf_complex convoluted_prs_time[len];
        fftwf_execute_dft(coarse_plan, &convoluted_prs[0], &convoluted_prs_time[0]);

        float tempVal;
        float maxVal=-99999;
        for (s=0;s<len;s++) {
            tempVal = sqrtf((convoluted_prs_time[s][0] * convoluted_prs_time[s][0]) + (convoluted_prs_time[s][1] * convoluted_prs_time[s][1]));
            if (tempVal>maxVal) {
                maxVal = tempVal;
            }
        }

        if (maxVal>global_max) {
            global_max = maxVal;
            global_max_pos = k;
        }
    }
    return global_max_pos;
}

double StandaloneEtiDecoder::get_fine_freq_corr(Csdr::complex<float> *input) {
    fftwf_complex *left;
    fftwf_complex *right;
    fftwf_complex *lr;
    double angle[504];
    double mean=0;
    double ffs;
    left = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * 504);
    right = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * 504);
    lr = (fftwf_complex*) fftwf_malloc(sizeof(fftwf_complex) * 504);
    uint32_t i;
    for (i = 0; i < 504; i++) {
        left[i][0] = input[2656 + 2048 + i].i();
        left[i][1] = input[2656 + 2048 + i].q();
        right[i][0] = input[2656 + i].i();
        right[i][1] = input[2656 + i].q();
    }
    for (i = 0; i < 504; i++){
        lr[i][0] = (left[i][0] * right[i][0] - left[i][1] * (-1)*right[i][1]);
        lr[i][1] = (left[i][0] * (-1)*right[i][1] + left[i][1] * right[i][0]);
    }

    for (i = 0; i < 504; i++){
        angle[i] = atan2f(lr[i][1],lr[i][0]);
    }
    for (i = 0; i < 504; i++){
        mean = mean + angle[i];
    }
    mean = (mean / 504);

    ffs = mean / (2 * M_PI) * 1000;

    fftwf_free(left);
    fftwf_free(right);
    fftwf_free(lr);

    return ffs;
}

void printUsage(const char* progname) {
    std::cerr << "Usage: " << progname << " [options]\n"
    << "Options:\n"
    << "  -i, --input <file>   Input file with IQ data (default: stdin)\n"
    << "  -o, --output <file>  Output file for ETI data (default: stdout)\n"
    << "  -h, --help           Show this help message\n";
}

int main(int argc, char* argv[]) {
    std::string input_filename;
    std::string output_filename;
    bool use_stdin = true;
    bool use_stdout = true;

    static struct option long_options[] = {
        {"input", required_argument, 0, 'i'},
        {"output", required_argument, 0, 'o'},
        {"help", no_argument, 0, 'h'},
        {0, 0, 0, 0}
    };

    int option_index = 0;
    int c;
    while ((c = getopt_long(argc, argv, "i:o:h", long_options, &option_index)) != -1) {
        switch (c) {
            case 'i':
                input_filename = optarg;
                use_stdin = false;
                break;
            case 'o':
                output_filename = optarg;
                use_stdout = false;
                break;
            case 'h':
                printUsage(argv[0]);
                return 0;
            default:
                printUsage(argv[0]);
                return 1;
        }
    }

    std::signal(SIGINT, handle_signal);
    std::signal(SIGTERM, handle_signal);


    // Open files
    std::ifstream input;
    if (!use_stdin) {
        input.open(input_filename, std::ios::binary);
        if (!input.is_open()) {
            std::cerr << "Error: Cannot open input file\n";
            return 1;
        }
    }
    std::istream& inputStream = use_stdin ? std::cin : input;

    std::ofstream output;
    if (!use_stdout) {
        output.open(output_filename, std::ios::binary);
        if (!output.is_open()) {
            std::cerr << "Error: Cannot open output file\n";
            return 1;
        }
    }
    std::ostream& outputStream = use_stdout ? std::cout : output;

    std::cerr << "Starting DAB ETI decoder (standalone)\n";

    try {
        StandaloneEtiDecoder decoder(outputStream);
        std::vector<float> iq_buffer(524288 * 2);
        size_t total_samples = 0;
        size_t frames_processed = 0;

        // while (inputStream.read(reinterpret_cast<char*>(iq_buffer.data()),
        // iq_buffer.size() * sizeof(float))) {
        while (!g_stop_requested &&
            inputStream.read(reinterpret_cast<char*>(iq_buffer.data()),
            iq_buffer.size() * sizeof(float))) {

        size_t bytes_read = inputStream.gcount();
        size_t samples_read = bytes_read / (2 * sizeof(float));

        if (samples_read == 0) break;

        total_samples += samples_read;

            // Convert to complex and add to decoder
            std::vector<Csdr::complex<float>> samples(samples_read);
            for (size_t i = 0; i < samples_read; ++i) {
                samples[i] = Csdr::complex<float>(
                    iq_buffer[i * 2],
                    iq_buffer[i * 2 + 1]
                );
            }

            decoder.addInputData(reinterpret_cast<uint8_t*>(samples.data()),
                                 samples.size() * sizeof(Csdr::complex<float>));

            // Process all available frames
            while (!g_stop_requested && decoder.hasEnoughData()) {
                decoder.processFrame();
                frames_processed++;
            }
        }

        //    std::cerr << "Decoding complete!\n";
        //    std::cerr << "Samples: " << total_samples << ", Frames: " << frames_processed << "\n";

            if (g_stop_requested) {
                std::cerr << "Interrupted, shutting down...\n";
            }
            std::cerr << "Decoding complete!\n";
            std::cerr << "Samples: " << total_samples << ", Frames: " << frames_processed << "\n";

    } catch (const std::exception& e) {
        std::cerr << "Error: " << e.what() << "\n";
        return 1;
    }

    return 0;
}
