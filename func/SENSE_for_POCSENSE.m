function im = SENSE_for_POCSENSE(kData, RecMtx)
    %Used for Recon SENSE = 1
    I = ifft2c(kData);
    im = sum(RecMtx.*I, 3);
end