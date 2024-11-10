import requests
import re
import shutil
from lxml import etree
from pathlib import Path


class DownloadCompose:
    table_link = "https://compose.obspm.fr/table"
    base_link = "https://compose.obspm.fr"
    target_files = ("eos.t", "eos.nb", "eos.thermo")

    def __init__(self, download_dir=Path("./downloads/compose")):
        self.download_dir = download_dir
        self.title_and_link = DownloadCompose.get_title_and_link()

    def get_title_and_link():
        print(f"DownloadCompose: Fetching data from {DownloadCompose.table_link} \n...")
        r = requests.get(DownloadCompose.table_link)
        # print(f"{r.status_code=}")
        # print(f"{r.encoding=}")
        tree = etree.HTML(r.text)
        odd_line = tree.xpath("//tr[@class='odd']")
        even_line = tree.xpath("//tr[@class='even']")
        all_line = odd_line + even_line
        print(
            f"DownloadCompose: Find {len(all_line)} EOS data on website {DownloadCompose.table_link}"
        )
        ret = {}
        for line in all_line:
            title = line.xpath("td[2]/text()")[0].strip()
            link = line.xpath(".//@href")[0]
            id = int(re.search(r"\d+", link).group())
            if id in ret:
                raise ValueError(f"Duplicate id={id}")
            ret[id] = (title, link)
        return dict(sorted(ret.items()))

    def print_eos_list(self):
        for id, (title, _) in self.title_and_link.items():
            print(f"id = {id:3}, name = {title}")

    def eos_name(self, id: int):
        return self.title_and_link[id][0]

    def eos_download_dir(self, id: int):
        dir = self.download_dir / f"{id:03}"
        dir.mkdir(parents=True, exist_ok=True)
        return dir

    def requests_download(url: str, folder, force=False):
        local_filename = Path(folder) / url.split("/")[-1]
        if local_filename.exists() and not force:
            print(f"{local_filename} already exists, skip download {url}")
            return

        with requests.get(url, stream=True, allow_redirects=True) as r:
            with open(local_filename, "wb") as f:
                shutil.copyfileobj(r.raw, f)
        print(f"Downloaded {url} to {local_filename}")

    def download_id(self, id: int, force=False):
        _, tail_link = self.title_and_link[id]
        eos_link = DownloadCompose.base_link + tail_link
        r = requests.get(eos_link)
        tree = etree.HTML(r.text)
        hrefs = tree.xpath("//a/@href")
        # print(f"{hrefs=}")
        target_hrefs = [
            href for href in hrefs if href.endswith(DownloadCompose.target_files)
        ]
        # print(f"{target_hrefs=}")
        for href in target_hrefs:
            url = DownloadCompose.base_link + href
            DownloadCompose.requests_download(
                url,
                self.eos_download_dir(id),
                force=force,
            )
